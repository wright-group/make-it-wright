import WrightTools as wt
import pathlib

from psutil import virtual_memory

from .andor import fromAndorNeo
from .gwyddion import fromGwyddion_traces
from .sp130 import fromSP130
from .horiba import fromLabramHR, horiba_typeID
from .iontof import fromITA, ITApeaks
from .xrd import fromBruker


px_per_um = {
    '5x-Jin' : 0.893,
    '20x-Jin' : 3.52, 
    '100x-Wright' : 18.2,
    '5' : 0.893,
    '20' : 3.52,
    '100' : 18.2,
}
# add 10x for now with approximation
px_per_um["10"] = px_per_um["10x-Jin"] = 2 * px_per_um["5"]


def typeID(*fpaths):
    """
    Infer what kind of data the file contains.
    The kind will inform on how to correctly import the data.
    """
    types = {}
    for fpath in map(pathlib.Path, fpaths):
        if fpath.suffix == '.ita':
            types[fpath] = 'iontof_SIMS'
            print(f"file {fpath} is IonToF SIMS data")
            
        if fpath.suffix == '.txt':
            with open(fpath) as f:
                txt = f.read()
            if "LabRAM HR" in txt:
                if (htype := horiba_typeID(fpath)) is not None:
                    types[fpath] = htype
            if "Goniometer" in txt:
                types[fpath] = 'Bruker_XRD'
            if "[m]" in txt:
                types[fpath] = 'Gwyddion_traces'

        if fpath.suffix == '.asc':
            with open(fpath) as f:
                txt = f.read()
            if "*BLOCK" in txt:
                types[fpath] = 'TRPL'
            else:
                types[fpath] = 'ASCII'

        if fpath.suffix == '.wt5':
            types[fpath] = 'wt5'

    print(f"{len(types)} of {len(fpaths)} files identified as valid data types")
    return types


def listfiles(fdir:str|pathlib.Path, pattern:str="*") -> list[pathlib.Path]:
    """Generate a list of filepaths within a directory.
    Includes files from nested directories.

    Parameters
    ----------
    fdir: path-like
        directory to walk
    pattern: string
        pattern used to filter files.  default uses no filter
    """
    return [
        pi for pi in filter(
            lambda pi: pi.is_file(), pathlib.Path(fdir).rglob(pattern)
        )
    ]
    

def parse(fdir, objective, select_types=None, keywords:list|str=[], exclude=[]):
    """
    DOCUMENTATION NEEDED
    """
    files = listfiles(fdir)
    
    include = [1 for i in range(len(files))]
    if keywords:
        if type(keywords) is not list:
            keywords = [keywords]
        for kw in keywords:
            for i, f in enumerate(files):
                if kw not in str(f):
                    include[i]=0
    if exclude:
        if type(exclude) is not list:
            exclude = [exclude]
        for x in exclude:
            for i, f in enumerate(files):
                if x in str(f):
                    include[i]=0

    files = [file for i, file in zip(include, files) if i]
    print(f'found {sum(include)} files matching keyword specifications')
    
    ftypes = typeID(*files)
    if select_types:
        to_delete=[]
        num_removed=0
        for key, value in ftypes.items():
            if value not in select_types:
                to_delete.append(key)
                num_removed+=1
        if to_delete:
            for key in to_delete:
                del(ftypes[key])
            print(f'excluded {num_removed} files that did not match specified data type(s)')
    
    if 'ASCII' in ftypes.values():
        if not objective:
            objective = input(f'Enter objective lens magnification if all data in this directory used the same lens. Otherwise, press enter: ')
        if not objective:
            objective = 'prompt'

    #make sure call doesn't generate too much data, roughly 1 GB
    too_much_data = False
    if len([dtype for dtype in ftypes.values() if dtype=='iontof_SIMS']) > 1:
        too_much_data = True
    if len([dtype for dtype in ftypes.values() if dtype=='ASCII']) > 100:
        too_much_data = True
    if len(ftypes) > 200:
        too_much_data = True
    if sum([f.stat().st_size for f in files]) > virtual_memory().available:
        too_much_data = True

    if too_much_data:
        raise MemoryError("too much data in directory, parsing cancelled to prevent storage overflow")

    d = []
    for fpath, dtype in ftypes.items():
        basename = fpath.stem

        if dtype.startswith('LabramHR'):
            d.append(fromLabramHR(fpath, name=basename))

        elif dtype=='Bruker_XRD':
            l0 = len(d)
            d = d + fromBruker(fpath)

        elif dtype=='Gwyddion_traces':
            d.append(fromGwyddion_traces(fpath, name=None, ID_steps=True))

        elif dtype=='iontof_SIMS':
            d.append((fpath, ITApeaks(fpath)))

        elif dtype=='TRPL':
            l0 = len(d)
            d.append(fromSP130(fpath, name=basename))
            print(basename)

        elif dtype=='ASCII':
            try:
                d.append(fromAndorNeo(fpath, name=basename, px_per_um=px_per_um[objective] if objective else None))
            except Exception as e:
                print(f'attempted to extract ASCII data from path <{fpath}> but it was not recognized by the andor module')
                raise e
            print(basename)
        
        elif dtype=='wt5':
            d.append(wt.open(fpath))
    if len(d)==1:
        d=d[0]
    return d
