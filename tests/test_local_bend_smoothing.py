import tempfile
import unittest
import json
from unittest.mock import patch
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from pipeline.local_bend_smoothing import (ATTRIBUTES, OUTPUT_ATTRIBUTES, prepare_bends,
    build_topology, candidate_weights, smooth_local)
from pipeline.spatial_smoothing import smooth_attributes, bend_neighbor_graph, run_spatial_smoothing
from pipeline.clustering_confinement import open_dataset_confinement_clustering


def chain(lengths):
    n=len(lengths)
    df=pd.DataFrame(dict(file=['oc']*n,networkGraph=[1]*n,combined_reach_id=[10]*n,
        combined_reach_up=[np.nan]*n,combined_reach_dn=[np.nan]*n,
        bendDistOut=np.arange(n,0,-1),bendLen=lengths))
    for i,a in enumerate(ATTRIBUTES):
        df[a]=np.arange(n,dtype=float)+i
    return prepare_bends(df)


class LocalSmoothingTests(unittest.TestCase):
    def test_equal_length_weights_and_candidates(self):
        df=chain([100]*7)
        topology=build_topology(df,4)
        ids,actual,effective,w,directions=candidate_weights(topology,3,3,.5,True)
        np.testing.assert_array_equal(ids,[3,2,1,0,4,5,6])
        np.testing.assert_allclose(actual,[0,100,200,300,100,200,300])
        np.testing.assert_array_equal(actual,effective)
        expected=np.exp(-2*np.array([0,1,2,3,1,2,3])**2)
        np.testing.assert_allclose(w,expected/expected.sum())

    def test_short_mirrored_neighbors(self):
        df=chain([100,100,25,100,25,100,100])
        topology=build_topology(df,3)
        _,actual,effective,w,_=candidate_weights(topology,3,3,.5,False)
        np.testing.assert_allclose(actual,[0,62.5,125,225,62.5,125,225])
        self.assertAlmostEqual(w[0],.4990968,places=5)
        _,actual,effective,w,_=candidate_weights(topology,3,3,1,True)
        np.testing.assert_allclose(effective,[0,100,200,300,100,200,300])
        self.assertAlmostEqual(w[0],.39905,places=4)

    def test_long_neighbor_and_boundary(self):
        df=chain([100,200,25])
        topology=build_topology(df,4)
        ids,actual,effective,_,directions=candidate_weights(topology,0,4,1,True)
        np.testing.assert_array_equal(ids,[0,1,2])
        np.testing.assert_allclose(actual,[0,150,262.5])
        np.testing.assert_allclose(effective,[0,150,300])
        self.assertEqual(directions,['self','down','down'])

    def test_cross_reach_and_network_boundary(self):
        df=chain([100]*4)
        df['combined_reach_id']=[10,10,20,20]
        df['combined_reach_dn']=[20,20,np.nan,np.nan]
        df['combined_reach_up']=[np.nan,np.nan,10,10]
        df=prepare_bends(df)
        topology=build_topology(df,3)
        np.testing.assert_array_equal(topology.paths[0,1],[1,2,3])
        df.loc[df.combined_reach_id==20,'networkGraph']=2
        topology=build_topology(prepare_bends(df),3)
        np.testing.assert_array_equal(topology.paths[0,1],[1,-1,-1])
        self.assertEqual(topology.diagnostics['missing_or_outside_network_references'],2)

    def test_original_arithmetic_parity(self):
        df=chain([100]*5)
        with tempfile.TemporaryDirectory() as tmp:
            ld=bend_neighbor_graph(df,Path(tmp)/'dist.pkl')
        expected=smooth_attributes('10_3',ATTRIBUTES,ld,df,max_neighbors=5)
        actual=smooth_local(df,neighbors=2,alpha=.5,length_floor=False)
        np.testing.assert_allclose(actual.loc[actual.bendID=='10_3',OUTPUT_ATTRIBUTES].to_numpy(dtype=float)[0],expected.astype(float),rtol=1e-13)

    def test_singleton_nan_invalid_lengths_and_cycles(self):
        df=chain([100])
        df.loc[0,ATTRIBUTES[0]]=np.nan
        output=smooth_local(df)
        self.assertTrue(np.isnan(output.loc[0,OUTPUT_ATTRIBUTES[0]]))
        np.testing.assert_array_equal(output[OUTPUT_ATTRIBUTES[4:]].values,np.zeros((1,4)))
        for length in [0,-1,np.nan,np.inf]:
            with self.assertRaises(ValueError):
                build_topology(chain([length]))
        df=chain([100,100])
        df['combined_reach_up']=10
        df['combined_reach_dn']=10
        with self.assertRaisesRegex(ValueError,'repeated'):
            build_topology(df)

    def test_nan_propagation_and_invalid_attributes(self):
        df=chain([100]*3)
        df.loc[0,ATTRIBUTES[0]]=np.nan
        output=smooth_local(df,neighbors=2)
        self.assertTrue(output[OUTPUT_ATTRIBUTES[0]].isna().all())
        self.assertTrue(output[OUTPUT_ATTRIBUTES[4]].isna().all())
        df.loc[0,ATTRIBUTES[0]]=np.inf
        with self.assertRaisesRegex(ValueError,'Infinite'):
            smooth_local(df)

    def test_order_independence_and_parameter_validation(self):
        df=chain([100]*5)
        shuffled=prepare_bends(df.sample(frac=1,random_state=1))
        pd.testing.assert_frame_equal(smooth_local(df),smooth_local(shuffled))
        topology=build_topology(df,4)
        for n,alpha in [(0,.5),(5,.5),(2.5,.5),(2,0),(2,np.nan)]:
            with self.assertRaises(ValueError):
                candidate_weights(topology,0,n,alpha)
        df.loc[0,'bendDistOut']=df.loc[1,'bendDistOut']
        with self.assertRaisesRegex(ValueError,'ambiguous'):
            prepare_bends(df)

    def test_experiment_isolation_and_step8_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            values=root/'single_values'
            values.mkdir()
            xr.Dataset({'placeholder':('index',[1])}).to_netcdf(values/'global_50_02_conf.nc')
            config=root/'paths.json'
            config.write_text(json.dumps({'results_root':str(root)}))
            out=root/'experiment'
            df=chain([100]*5)
            with patch('pipeline.spatial_smoothing._prepare_smoothing_dataframe',return_value=df):
                result=run_spatial_smoothing(config_path=config,workers=1,method='local',output_dir=out)
            loaded,path=open_dataset_confinement_clustering(2,config_path=config,input_dir=out)
            self.assertEqual(len(loaded),5)
            self.assertEqual(path,result['global_output'])
            with self.assertRaisesRegex(ValueError,'different smoothing'):
                run_spatial_smoothing(config_path=config,method='local',output_dir=out,alpha=1)
            self.assertFalse((root/'single_smoothed').exists())

    def test_canonical_local_directory_rejects_legacy_and_unmarked_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            values=root/'single_values'
            values.mkdir()
            xr.Dataset({'placeholder':('index',[1])}).to_netcdf(values/'global_50_02_conf.nc')
            config=root/'paths.json'
            config.write_text(json.dumps({'results_root':str(root)}))
            destination=root/'single_smoothed'
            destination.mkdir()
            (destination/'length_dict_af.pkl').write_bytes(b'incomplete legacy output')
            with self.assertRaisesRegex(ValueError,'without a settings manifest'):
                run_spatial_smoothing(config_path=config,workers=1,method='local')
            (destination/'length_dict_af.pkl').unlink()
            with patch('pipeline.spatial_smoothing._prepare_smoothing_dataframe',return_value=chain([100]*5)):
                result=run_spatial_smoothing(config_path=config,workers=1)
            self.assertEqual(result['global_output'],(destination/'global_50_02_smoothed.nc').resolve())
            with self.assertRaisesRegex(ValueError,'different smoothing'):
                run_spatial_smoothing(config_path=config,workers=1,method='legacy')

    def test_canonical_global_rejects_partial_continent_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            values=root/'single_values'
            values.mkdir()
            xr.Dataset({'placeholder':('index',[1])}).to_netcdf(values/'global_50_02_conf.nc')
            config=root/'paths.json'
            config.write_text(json.dumps({'results_root':str(root)}))
            df=chain([100]*5)
            df.loc[0,'file']='af'
            with patch('pipeline.spatial_smoothing._prepare_smoothing_dataframe',return_value=df):
                with self.assertRaisesRegex(ValueError,'requires all available continents'):
                    run_spatial_smoothing(config_path=config,workers=1,continents=['oc'])


if __name__=='__main__':
    unittest.main()
