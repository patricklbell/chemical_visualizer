import { ReactComponent as ArrowForwardIcon } from 'assets/ArrowForward.svg';

import { useMolInfo } from "hooks/useWasm";
import { Accordian, AccordianButton, AccordianPanel } from "components/Accordian";

const MolStructure = ({viewer, refresh}) => {
    const [info] = useMolInfo(viewer, refresh);

    if (!info)
      return <></>

    return <div className="flex flex-col gap-4">
      {info?.title.replace( /\s/g, '') && <h1>{info?.title}</h1>}
      {info?.comments.replace( /\s/g, '') && <h1>{info?.comments}</h1>}
      {info?.metadata.replace( /\s/g, '') && <h1>{info?.metadata}</h1>}

      {info?.atoms.length ?
        <div className="border-2 border-accent1 p-1 mt-[0.6rem] pt-[0.6rem] rounded-xl relative">
          <span className="bg-back absolute left-3 top-[-1.0rem] px-1 text-accent font-semibold">Atoms</span>
          
          <Accordian>
            <AccordianButton>Atoms</AccordianButton>

            <AccordianPanel as="ul">
              <li key={-1} className="flex flex-row gap-5 justify-between border-b border-back1 mb-2">
                  <span className="w-[5%]">Symbol</span>
                  <span className="w-1/4">Charge</span>
                  <span className="w-1/4">Valence</span>
              </li>
              {info?.atoms.map((r, i) => 
              <li key={i} className="flex flex-row gap-5 justify-between">
                  <span className="w-[5%]">{r?.symbol}</span>
                  <span className="w-1/4">{r?.charge}</span>
                  <span className="w-1/4">{r?.valence}</span>
              </li>)}
            </AccordianPanel>
          </Accordian>
        </div> : <></>}

      {info?.bonds.length ? 
        <div className="border-2 border-bccent1 p-1 mt-[0.6rem] pt-[0.6rem] rounded-xl relative">
          <span className="bg-back absolute left-3 top-[-1.0rem] px-1 text-bccent font-semibold">Bonds</span>
          <Accordian>
            <AccordianButton>Bonds</AccordianButton>

            <AccordianPanel as="ul">
              <li key={-1} className="flex flex-row gap-5 justify-between border-b border-back1 mb-2">
                  <span className="w-1/4">Atom 1</span>
                  <ArrowForwardIcon className="h-4"/>
                  <span className="w-1/4 text-right pr-1">Atom 2</span>
              </li>
              {info?.bonds.map((b, i) => 
              <li key={i} className="flex flex-row gap-5 justify-between">
                  <span className="w-1/4">{b?.atom1Id}</span>
                  <ArrowForwardIcon className="h-4"/>
                  <span className="w-1/4 text-right pr-1">{b?.atom2Id}</span>
              </li>)}
            </AccordianPanel>
          </Accordian>
        </div> : <></>}
    </div>
}

export default MolStructure;