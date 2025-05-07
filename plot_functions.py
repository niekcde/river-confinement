import matplotlib.pyplot as plt
import numpy as np

from calc_functions import confinement_values

def plot_confinemend_apex(crid, po, pi,cdo, cdi, bendWidth, confiningFactor = ((1/28)*2) ):
    (poIntercept, piIntercept,
     centerPointHeight, confinementHeight,
     slopeOut, slopeInn, EROut, ERInn) = confinement_values(po, pi, cdo, cdi, 
                                                            bendWidth, confiningFactor)


    minH = min(min(pi), min(po))
    maxH = max(max(pi), max(po))


    #####
    # plot settings:
    axis_label_size = 20
    axis_tick_size  = 16

    outsideColor = 'silver'
    innsideColor = 'dimgray'


    ############
    # Figure 1
    ############
    f, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = [10,8])
    # f.suptitle(ID)
    ax1.plot(cdo             , po, color = 'black')
    ax1.plot(np.array(cdi)*-1, pi, color = 'black')

    ##################### 
    # Slope lines
    ax1.plot([0, poIntercept], [centerPointHeight, confinementHeight], color = 'blue')
    ax1.plot([0, piIntercept*-1], [centerPointHeight, confinementHeight], color = 'blue')

    #####################
    # River extend lines
    ax1.vlines(bendWidth/2  , minH, maxH, linestyle = '--', color = 'blue')
    ax1.vlines((bendWidth/2)*-1, minH, maxH, linestyle = '--', color = 'blue')


    #####################
    # boundary height line
    ax1.hlines(confinementHeight, max(cdo)*-1, max(cdi), linestyle = '--', color = 'red')

    #####################
    # River centerpoint
    ax1.scatter(0, centerPointHeight, marker = '*', c = 'red', s = 30
                , zorder = 100)

    loc = [0.6,1]
    
    ax1.text(cdo[-1]*loc[0] , maxH *loc[1], "Outside", fontsize=axis_label_size, color=outsideColor, ha="center")
    ax1.text(cdi[-1]*-loc[0], maxH *loc[1], "Innside", fontsize=axis_label_size, color=innsideColor, ha="center")

    ax1.set_xlabel('Distance (m)'  , fontsize = axis_label_size)
    ax1.set_ylabel('Elevation (m)' , fontsize = axis_label_size)

    # Customize tick size and color
    ax1.tick_params(axis='both', which='major', labelsize=axis_tick_size, colors='black')  # Major ticks

    plt.show()

    ###########
    # Figure 2
    ###########
    # f, ax2 = plt.subplots(nrows = 1, ncols = 1, figsize = [10,8])
    # f.suptitle(ID)
    # pointSize = 150


    # ax2.plot(*combinedLine.xy, zorder =10, color = 'red', linewidth = 3, label = 'Reach')
    # plotConnection(dfReach,reachCRS, df, ax2)

    # for i, ap in enumerate(apexP):
    #     if i != apexID:
    #         ax2.scatter(*ap.xy, c = 'black', marker = '^', zorder = 100, s= pointSize)
    # ax2.scatter(*apexP.xy, c = 'black', marker = '*', zorder = 1000, s= pointSize)

    # ax2.plot(*lineOut.xy, color = outsideColor, label = 'outside')
    # ax2.plot(*lineInn.xy, color = innsideColor, label = 'innside')



    # ax2.set_aspect('equal')
    # ax2.axis('off')

    # # plot Raster
    # rasterValues = raster[0,:,:]
    # rasterValues = rasterValues.values.ravel()
    # vmin = rasterValues[rasterValues != -9999].min()
    # im = raster.plot(ax = ax2, vmin = vmin, vmax = raster.max(), cmap = 'terrain', zorder = 0,
    #             add_colorbar=False)
    # ax2.set_title('')



    # divider = make_axes_locatable(ax2)
    # cax = divider.append_axes("right", size="5%", pad=0.05)

    # cbar = plt.colorbar(im, cax=cax, shrink=0.8, aspect=10)

    # # Set the label with a custom text size
    # # cbar.set_label('Colorbar Label', fontsize=14)  # Adjust label text size
    # # Adjust the size of the tick labels
    # cbar.ax.tick_params(labelsize=axis_tick_size)  # Adjust tick label text size
    # # plt.colorbar(im, cax=cax, fontsize = axis_tick_size)

    # ax2.legend()


    # plt.tight_layout()
    # plt.savefig(directory + 
    #             f'results/figures/methods_1_{c}_{i}_{dfReach.iloc[0].combined_reach_id}.png')
    # plt.show()