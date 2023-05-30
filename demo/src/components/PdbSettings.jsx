import Combobox from 'components/Combobox';
import PrimaryButton from 'components/PrimaryButton';
import Check from 'components/Check';

import { usePdbDrawSettings } from "../hooks/useWasm";
import { useState } from 'react';

const PdbSettings = ({ viewer, setLoaded = () => {}, setLoadingMessage = () => {} }) => {
    const [pdbDrawSettings, setPdbDrawSettings] = usePdbDrawSettings(viewer);
    const [dictionaryLoaded, setDictionaryLoaded] = useState(false);

    return <>
        <span className="text-center text-[1.2rem] font-semibold">PDB Settings</span>
        {(() => {
            const options = ["Polypeptide Chains", "Secondary Structures", "Amino Acids"]
            const modes = {
                [viewer.ResidueColorMode.CHAIN.value]: "Polypeptide Chains",
                [viewer.ResidueColorMode.SECONDARY.value]: "Secondary Structures",
                [viewer.ResidueColorMode.AMINO_ACID.value]: "Amino Acids",
                "Polypeptide Chains": viewer.ResidueColorMode.CHAIN,
                "Secondary Structures": viewer.ResidueColorMode.SECONDARY,
                "Amino Acids": viewer.ResidueColorMode.AMINO_ACID,
            };

            return <div className='flex flex-row items-center'>
                <span className='pr-4'>Residue Colouring</span>
                <Combobox 
                    className="flex-grow"
                    selected={modes[pdbDrawSettings?.residueColorMode.value]}
                    setSelected={(opt) => 
                        setPdbDrawSettings({
                        ...pdbDrawSettings,
                        residueColorMode: modes[opt],
                        })
                    }
                    options={options}
                />
            </div>
        }
        )()}

        <div className="flex flex-row items-center">
            <Check
                value={pdbDrawSettings?.drawHeteroAtoms}
                setValue={(v) => setPdbDrawSettings({...pdbDrawSettings, drawHeteroAtoms: v})}
            />
            <span className="px-4">Show Heterogen Atoms</span>
            </div>

            <div className="flex flex-row items-center">
            <Check
                value={pdbDrawSettings?.drawResidueAtoms}
                setValue={(v) => {
                    if (!dictionaryLoaded && v) {
                        setLoaded(false);
                        setLoadingMessage("Loading Residue Dictionary");
                        fetch('examples/het_dictionary.pdb')
                            .then(res => res.text())
                            .then(txt => {
                                viewer.loadPdbDictionary(txt);
                                setDictionaryLoaded(true);
                                setLoaded(true);
                            });
                        setPdbDrawSettings({...pdbDrawSettings, drawResidueAtoms: v});
                    } else {
                        setPdbDrawSettings({...pdbDrawSettings, drawResidueAtoms: v});
                    }
                }}
            />
            <span className="px-4">Show Atoms Inside Chains</span>
            </div>

            <div className="flex flex-row items-center">
            <Check
                value={pdbDrawSettings?.drawResidueRibbons}
                setValue={(v) => setPdbDrawSettings({...pdbDrawSettings, drawResidueRibbons: v})}
            />
            <span className="px-4">Show Ribbons</span>
            </div>

            
            <div className="flex flex-row items-center">
            <Check
                value={pdbDrawSettings?.drawWaterAtoms}
                setValue={(v) => setPdbDrawSettings({...pdbDrawSettings, drawWaterAtoms: v})}
            />
            <span className="px-4">Show Water Atoms</span>  
        </div>
    </>
}

export default PdbSettings;