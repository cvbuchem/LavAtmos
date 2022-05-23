from matplotlib import pyplot as plt
import matplotlib.patheffects as pe

fig_dir = '/home/jovyan/ThermoEngine/LavAtmos/figs/'

def example1_plot(T, lavatmos_bse):

    fig, axes = plt.subplots(1,2,figsize=(21,14))

    # Define included elements
    elements = ['O', 'Si','Al','Ti','Fe','Mg','Ca','Na','K']
    
    # Define color scheme
    colors_mokole1 = ['blue','dodgerblue','aqua','green','lime','gold','red','rosybrown','deeppink']
    color_dict_mokole = {}
    for i,el in enumerate(elements):
        color_dict_mokole[el] = colors_mokole1[i]
        
    # Determine linestyles
    linestyles = ['solid','dashed','dashdot','dotted',(0, (3, 1, 1, 1, 1, 1)),(0,(1,10))]
    linestyle_dict = {}
    linewidth = 3

    # Plot parameters
    axis_labelsize = 20
    tick_label_fontsize = 16
    title_size = 20

    counter = {}

    for spec in lavatmos_bse.columns:
        # if spec not in included_species:
            # print(spec)
            # continue
        # print(spec)
        if spec[:1] in elements:
            el = spec[:1]
        elif spec[:2] in elements:
            el = spec[:2]
        try:
            counter[el] += 1
        except:
            counter[el] = 0

        color = color_dict_mokole[el]
        linestyle_dict[spec] = linestyles[counter[el]] 
        linestyle = linestyle_dict[spec]
        
        # Partial pressure
        axes[0].plot(T,lavatmos_bse[spec],label=spec,\
                     color=color,linestyle=linestyle,linewidth=linewidth,\
                     path_effects=[pe.Stroke(linewidth=linewidth+1, foreground='darkgrey'), pe.Normal()])
        
        # Mole fractions
        axes[1].plot(T,lavatmos_bse[spec]/lavatmos_bse.sum(axis=1),label=spec,\
                     color=color,linestyle=linestyle,linewidth=linewidth,\
                     path_effects=[pe.Stroke(linewidth=linewidth+1, foreground='darkgrey'), pe.Normal()])
        
    # Settings for partial pressure plot
    axes[0].set_title('Partial pressure',fontsize=title_size)
    axes[0].set_ylabel('Partial pressure [bar]',fontsize=axis_labelsize)
    axes[0].set_xlabel('Temperature [K]',fontsize=axis_labelsize)
    axes[0].set_ylim(1e-10,100)
    axes[0].set_xlim(1500,4000)
    axes[0].set_yscale('log')
    axes[0].grid(visible=True, which='major', axis='both', zorder=0, alpha=0.3)
    axes[0].tick_params(axis='both',which='both',labelsize=tick_label_fontsize)

    # Settings for mole faction plot
    axes[1].set_title('Molefraction',fontsize=title_size)
    axes[1].set_ylabel('Mole fraction',fontsize=axis_labelsize)
    axes[1].set_xlabel('Temperature [K]',fontsize=axis_labelsize)
    axes[1].set_ylim(1e-7,1)
    axes[1].set_xlim(1500,4000)
    axes[1].set_yscale('log')
    axes[1].grid(visible=True, which='major', axis='both', zorder=0, alpha=0.3)
    axes[1].tick_params(axis='both',which='both',labelsize=tick_label_fontsize)

    axes[1].legend(bbox_to_anchor=(1.05, 0.95), loc='upper left', fontsize=16, frameon=False)

    plt.subplots_adjust(wspace=0.2, hspace=1)
    plt.savefig(fig_dir+'example1_results_BSE.jpg',bbox_inches="tight")

    plt.show()