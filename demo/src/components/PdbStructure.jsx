import { usePdbInfo } from "hooks/useWasm"
import { Accordian, AccordianButton, AccordianPanel } from "components/Accordian";

const PdbStructure = ({viewer, refresh}) => {
    const [pdbInfo] = usePdbInfo(viewer, refresh)

    return <ul>
    {pdbInfo?.models.map((m, i) => 
      <Accordian as="li" key={i} defaultOpen={m?.id === '0'}>
        <AccordianButton>
          Model {m?.id}
        </AccordianButton>

        <AccordianPanel className="pl-4">
          {m?.chains.length ? 
            <div className="border-2 border-accent1 p-1 mt-[0.6rem] pt-[0.6rem] rounded-xl relative">
              <span className="bg-back absolute left-3 top-[-1.0rem] px-1 text-accent font-semibold">Primary Structures</span>
              
              <Accordian>
                <AccordianButton>Polypeptide Chains</AccordianButton>

                <AccordianPanel as="ul">
                  {m.chains.map((c, i) => 
                  <Accordian as="li" key={i}>
                    <AccordianButton>Chain {c?.id}</AccordianButton>

                    <AccordianPanel as="ul">
                      <li key={-1} className="flex flex-row gap-5 justify-between border-b border-back1 mb-2">
                        <span className="w-[5%]">Id</span>
                        <span className="w-1/2">Name</span>
                      </li>
                      {c?.residues.map((r, i) => 
                      <li key={i} className="flex flex-row gap-5 justify-between">
                        <span className="w-[5%]">{r?.id}</span>
                        <span className="w-1/2">{r?.name}</span>
                      </li>)}
                    </AccordianPanel>
                  </Accordian>)}
                </AccordianPanel>
              </Accordian>
            </div> : <></>
          }
        </AccordianPanel>

        <AccordianPanel className="pl-4">
          {m?.sheets.length || m?.helices.length ? 
            <div className="border-2 border-bccent1 p-1 mt-[0.6rem] pt-[0.6rem] rounded-xl relative">
              <span className="bg-back absolute left-3 top-[-1.0rem] px-1 text-bccent font-semibold">Secondary Structures</span>
              {m?.helices.length ?
              <Accordian>
                <AccordianButton>Helices</AccordianButton>

                <AccordianPanel as="ul">
                  <li key={-1} className="flex flex-row gap-5 justify-between items-center border-b border-back1 mb-2">
                    <span className="w-[5%]">Id</span>
                    <span className="w-1/5">Name</span>
                    <span className="w-1/2">Length (Residues)</span>
                  </li>
                  {m.helices.map((h, i) => 
                  <>
                    <li key={i} className="flex flex-row gap-5 justify-between items-center">
                      <span className="w-[5%]">{h?.id}</span>
                      <span className="w-1/5">{h?.name.replace( /\s/g, '')}</span>
                      <span className="w-1/2">{h?.length}</span>
                    </li>
                    {h?.comment.replace( /\s/g, '').length ? 
                      <li className="pl-5">{h?.comment}</li> : <></>
                    }
                  </>)}
                </AccordianPanel>
              </Accordian> : <></>}

              {m?.sheets.length ?
              <Accordian>
                <AccordianButton>Sheets</AccordianButton>

                <AccordianPanel as="ul">
                  {m.sheets.map((s, i) => 
                  <Accordian as="li" key={i}>
                    <AccordianButton>Sheet {s?.id}</AccordianButton>

                    <AccordianPanel as="ul">
                      <li key={-1} className="flex flex-row gap-5 justify-between border-b border-back1 mb-2">
                        <span className="w-[5%]">Id</span>
                        <span className="w-1/2">Sense</span>
                      </li>
                      {s?.strands.map((s, i) => 
                      <li key={i} className="flex flex-row gap-5 justify-between">
                        <span className="w-[5%]">{s?.id}</span>
                        <span className="w-1/2">{s?.sense === 0 && "Start"}{s?.sense === 1 && "Parallel"}{s?.sense === -1 && "Anti-parallel"}</span>
                      </li>)}
                    </AccordianPanel>
                  </Accordian>)}
                </AccordianPanel>
              </Accordian> : <></>}
            </div> : <></>
          }
        </AccordianPanel>
      </Accordian>
    )}
  </ul>
}

export default PdbStructure;