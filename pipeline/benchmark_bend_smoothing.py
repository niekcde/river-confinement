"""Benchmark each smoothing method in a fresh process on identical prepared input."""
if __package__ in (None, ''):
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    __package__='pipeline'

import argparse
import json
import multiprocessing as mp
import platform
import resource
import time
from pathlib import Path

import pandas as pd
import psutil

from .build_smoothing_sensitivity import load_input
from .local_bend_smoothing import build_topology, smooth_local
from .spatial_smoothing import run_bend_smoothing


def worker(method, input_pickle, output, neighbors, alpha, length_floor):
    df=pd.read_pickle(input_pickle)
    start=time.perf_counter()
    if method=='legacy':
        for continent, frame in df.groupby('file'):
            run_bend_smoothing(continent,frame,single_smoothed_dir=output,cross_token='benchmark',hf_token='baseline')
    else:
        for continent, frame in df.groupby('file'):
            frame=smooth_local(frame,build_topology(frame,neighbors),neighbors=neighbors,alpha=alpha,length_floor=length_floor)
            frame.to_xarray().to_netcdf(output/f'{continent}_benchmark_local.nc')
    elapsed=time.perf_counter()-start
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes=peak if platform.system()=='Darwin' else peak*1024
    (output/'measurement.json').write_text(json.dumps(dict(method=method,bends=len(df),seconds=elapsed,
        peak_process_rss_bytes=peak_bytes,includes='graph construction, averaging, output writing; excludes shared input preparation',
        neighbors=neighbors,alpha=alpha,length_floor=length_floor),indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--continents',nargs='+',required=True)
    p.add_argument('--neighbors',type=int,default=3)
    p.add_argument('--alpha',type=float,default=.75)
    p.add_argument('--length-floor',action=argparse.BooleanOptionalAction,default=True)
    p.add_argument('--legacy-timeout',type=float,default=1800,help='Seconds before terminating baseline; not treated as completed timing')
    p.add_argument('--legacy-memory-gb',type=float,default=4,help='Stop baseline at this process RSS; recorded as incomplete')
    args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter()
    df=load_input(args.input,args.continents)
    df.groupby(['file','networkGraph']).size().rename('bends').to_csv(args.output_dir/'network_sizes.csv')
    prepared=args.output_dir/'prepared.pkl'
    df.to_pickle(prepared)
    (args.output_dir/'input_preparation.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,input=str(args.input)),indent=2))
    del df
    records=[]
    for method in ['local','legacy']:
        output=args.output_dir/method
        output.mkdir()
        process=mp.get_context('spawn').Process(target=worker,args=(method,prepared,output,args.neighbors,args.alpha,args.length_floor))
        process.start()
        reason = None
        started = time.perf_counter()
        observed_peak = 0
        while process.is_alive():
            process.join(.2)
            try:
                observed_peak = max(observed_peak, psutil.Process(process.pid).memory_info().rss)
            except psutil.NoSuchProcess:
                pass
            if method == 'legacy' and process.is_alive():
                if time.perf_counter()-started > args.legacy_timeout:
                    reason = 'timeout'
                    break
                if observed_peak > args.legacy_memory_gb * 1024**3:
                    reason = 'memory_limit'
                    break
        if process.is_alive():
            process.terminate()
            process.join()
            records.append(dict(method=method,status=reason,timeout_seconds=args.legacy_timeout,
                                memory_limit_gb=args.legacy_memory_gb, observed_peak_rss_bytes=observed_peak))
        elif process.exitcode:
            records.append(dict(method=method,status='failed',exitcode=process.exitcode))
        else:
            records.append(dict(status='complete',**json.loads((output/'measurement.json').read_text())))
        (args.output_dir/'benchmark.json').write_text(json.dumps(records,indent=2))
        print(records[-1],flush=True)
    prepared.unlink()


if __name__=='__main__':
    main()
