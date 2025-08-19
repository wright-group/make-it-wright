try:
    import pySPM
except ImportError:
    pass
import numpy as np
import WrightTools as wt


def open_ita(fpath):
    try:
        ita = pySPM.ITA(fpath)
    except ModuleNotFoundError:
        print("""
            ionTOF support is optional and was not specified at install.
            to work with iontof data, please install the optional dependencies
            `pip install git+https://github.com/wright-group/makeitwright.git[iontof]`
            """
        )
    return ita


def fromITA(fpath, name=None, select_channels=None):

    ita = open_ita(fpath)
    ita.show_summary()
    summ = ita.get_summary()
    
    xarr = np.linspace(0, summ['fov']*1e6, num=summ['pixels']['x'])
    yarr = np.linspace(0, summ['fov']*1e6, num=summ['pixels']['y'])
    scarr = np.linspace(1, int(summ['Scans']), num=int(summ['Scans']))
    charrs = {}
    if select_channels is not None:
        idxs = []
        for peak in summ['peaks']:
            if peak['id'] in select_channels or peak['assign'] in select_channels:
                idxs = idxs + [peak['id']]
        for idx in idxs:
            if summ['peaks'][idx]['assign']:
                chname = summ['peaks'][idx]['assign']
            elif summ['peaks'][idx]['desc']:
                chname = summ['peaks'][idx]['desc']
            else:
                chname = str(int(summ['peaks'][idx]['cmass'])) + 'mz'
            charr = ita.getImage(idx,0)        
            for i in range(1,len(scarr)):
                j = ita.getImage(idx,i)
                charr = np.dstack((charr,j))
            charrs[chname] = charr
            print("channel <" + chname + "> found")
    else:
        for peak in summ['peaks']:
            if peak['assign']:
                chname = peak['assign']
            elif peak['desc']:
                chname = peak['desc']
            else:
                chname = str(int(peak['cmass'])) + 'mz'
            idx = peak['id']
            charr = ita.getImage(idx,0)        
            for i in range(1,len(scarr)):
                j = ita.getImage(idx,i)
                charr = np.dstack((charr,j))
            charrs[chname] = charr
            print("channel <" + chname + "> found")
        
    d = wt.Data()
    d.create_variable(name='x', values=xarr[:,None,None], units='um')
    d.create_variable(name='y', values=yarr[None,:,None], units='um')
    d.create_variable(name='scan', values=scarr[None,None,:], units='s')
    for chname, charr in charrs.items():
        d.create_channel(name=chname, values=charr)
    d.transform('x','y','scan')

    return d


def ITApeaks(fpath):
    ita = open_ita(fpath)
    summ = ita.get_summary()
    return summ['peaks']
