import { useEffect, useRef, useState } from "react";
import { Panel, PanelGroup, PanelResizeHandle } from "react-resizable-panels";
import { useResizeDetector } from 'react-resize-detector';
import classnames from 'classnames';

import { ReactComponent as ChevronRightIcon } from 'assets/ChevronRight.svg';
import { ReactComponent as SettingsIcon } from 'assets/Settings.svg';
import { ReactComponent as HelpIcon } from 'assets/Help.svg';
import { ReactComponent as DragClickIcon } from 'assets/DragClick.svg';
import { ReactComponent as ZoomIcon } from 'assets/Zoom.svg';
import { ReactComponent as DarkModeIcon } from 'assets/DarkMode.svg';
import { ReactComponent as LightModeIcon } from 'assets/LightMode.svg';
import { ReactComponent as MoreVertIcon } from 'assets/MoreVert.svg';

import Throbber from "components/Throbber";
import PrimaryButton from "components/PrimaryButton";
import Combobox from "components/Combobox";
import Modal from "components/Modal";
import {Tooltip, TooltipContent, TooltipTrigger} from "components/Tooltip";
import PdbSettings from "components/PdbSettings";
import PdbStructure from "components/PdbStructure";
import MolStructure from "components/MolStructure";

import { useWasm, useCameraType, useDrawSettings } from "hooks/useWasm";
import { vec4 } from "utils/glmVector";

const RightPanel = ({viewer, isProtein, refresh, viewerLoaded, setViewerLoaded, setLoadingMessage}) => {
  return <div className={classnames("flex flex-col max-h-[100vh] flex-grow", {'z-10 opacity-10 pointer-events-none': !viewerLoaded})}>
    <div className="py-2 px-5 flex flex-col flex-grow gap-2 overflow-y-auto">
      <span className="text-center text-[1.2rem] font-semibold">Structure</span>
      {isProtein ? <PdbStructure viewer={viewer} refresh={refresh} /> : <MolStructure viewer={viewer} refresh={refresh} />}
    </div>

    {isProtein &&
    <div className="py-4 pt-2 px-5 flex flex-col gap-3 border-t border-back1">
      <PdbSettings viewer={viewer} setLoaded={setViewerLoaded} setLoadingMessage={setLoadingMessage} />
    </div>}
  </div>
}

const App = () => {
  const { width: canvasWidth, height: canvasHeight, ref: canvasRef } = useResizeDetector();
  const rightPanel = useRef(null);
  const fileUploadRef = useRef(null);
  const [isHelpOpen, setIsHelpOpen] = useState(false);
  const [isSettingsOpen, setIsSettingsOpen] = useState(false);
  const [rightCollapsed, setRightCollapsed] = useState(false);
  
  const examples = {
    "Human Deoxyhaemoglobin": "examples/4hhb.pdb",
    "Insulin Replacement": "examples/1bzv.pdb",
    "Micrococcus Lysodeikticus Catalase": "examples/1gwe.pdb",
    "Immunoglobulin": "examples/1igt.pdb",
    "E. Coli Beta-barrel Protein": "examples/2jmm.pdb",
    "Human Haemoglobin": "examples/4n7n.pdb",
    "Covid-19 Protease": "examples/6w63.pdb",
    "Fullerene": "examples/fullerene.mol",
    "Acetylene": "examples/acetylene.mol",
    "ATP": "examples/ATP.mol",
    "Caffeine": "examples/caffeine.mol",
    "Capsaicin": "examples/capsaicin.mol",
    "Chloroquine": "examples/chloroquine.mol",
    "Ethanol": "examples/ethanol.mol",
  }
  const [isProtein, setIsProtein] = useState(false);
  const [fullRightPanel, setFullRightPanel] = useState(false);
  const [viewerLoaded, setViewerLoaded] = useState(false);
  const [loadingMessage, setLoadingMessage] = useState("Loading Viewer");
  const [selected, setSelected] = useState(Object.keys(examples)[0]);
  const [openSelection, setOpenSelection] = useState(selected);

  const [refresh, setRefresh] = useState(false);
  const viewer = useWasm(canvasRef);
  const [cameraType, setCameraType] = useCameraType(viewer);
  const [drawSettings, setDrawSettings] = useDrawSettings(viewer);

  const updateViewerBackground = () => {
    if (!drawSettings)
      return;
    // @todo @hardcoded
    if (
      localStorage.getItem('color-theme') === 'light' ||
      !document.documentElement.classList.contains('dark')
    ) {
      setDrawSettings({
        ...drawSettings, 
        clearColor: vec4(drawSettings.clearColor, 1, 1, 1, 1)
      });
    } else {
      setDrawSettings({
        ...drawSettings, 
        clearColor: vec4(drawSettings.clearColor, 20/255, 22/255, 27/255, 1)
      });
    }
  }
  const toggleDarkMode = () => {
    if (
        localStorage.getItem('color-theme') === 'light' ||
        !document.documentElement.classList.contains('dark')
      ) {
        document.documentElement.classList.add('dark');
        localStorage.setItem('color-theme', 'dark');
    } else {
        document.documentElement.classList.remove('dark');
        localStorage.setItem('color-theme', 'light');
    }

    updateViewerBackground();
  }
  useEffect(() => {
    if (!viewer || !drawSettings || viewerLoaded)
      return;

    updateViewerBackground()
    setViewerLoaded(true);
  }, [viewer, drawSettings]);

  const togglePanel = (ref) => {
    const panel = ref.current;
    if (!panel) return;

    if (rightCollapsed) {
      panel.expand();
    } else {
      panel.collapse();
    }
  };

  const viewFile = (txt, path) => {
    if (path.endsWith('.pdb')) {
      viewer.loadPdb(txt);
      setIsProtein(true);
    } else {
      viewer.loadMol(txt);
      setIsProtein(false);
    }
    setViewerLoaded(true);
    setRefresh(!refresh);
  }

  useEffect(() => {
    setOpenSelection(selected);

    if (!viewer) return;

    const path = examples[selected];
    setViewerLoaded(false);
    setLoadingMessage("Loading " + selected)
    fetch(path)
      .then(res => res.text())
      .then(txt => viewFile(txt, path));
  }, [selected, viewer])

  useEffect(() => {
    viewer?.setViewportSize(canvasWidth, canvasHeight);
  }, [viewer, canvasHeight, canvasWidth]);

  return (
    <PanelGroup 
      className="bg-back text-fore shadow-fore1 fill-fore absolute left-0 top-0 w-[100vw] h-[100vh]" 
      autoSaveId="viewer"
      direction="horizontal"
    >
      <Panel defaultSize={75} minSize={50} className="flex flex-col border-r border-back1">
        <div className="flex flex-row justify-end pb-1 border-back1 border-b items-center content-center px-1">
          <Tooltip padding={5} placement={"right"}>
            <TooltipTrigger className="hidden lg:block my-auto pl-5" onClick={toggleDarkMode}>
              <LightModeIcon className="w-6 h-6 block dark:hidden"/>
              <DarkModeIcon className="w-6 h-6 hidden dark:block"/>
            </TooltipTrigger>
            
            <TooltipContent>
              <div className="bg-back p-2 border border-back1 rounded-md">
                <span className="block m-auto text-sm">
                  <span className="block dark:hidden text-fore">Light Mode</span>
                  <span className="hidden dark:block text-fore">Dark Mode</span>
                </span>
              </div>
            </TooltipContent>
          </Tooltip>

          <div className="flex-grow flex flex-row justify-center items-stretch gap-5">
            <span className="hidden lg:inline-flex items-center">Choose Example</span>
            <Combobox 
              className="text-sm lg:text-2xl text-center lg:w-[20rem]"
              inputClassName="font-bold"
              selected={openSelection}
              setSelected={setSelected}
              options={Object.keys(examples)}
            />
            <span className="inline-flex items-center text-fore1 text-opacity-50">Or</span>
            <PrimaryButton 
              onClick={() => {fileUploadRef.current.click()}}
              className="py-1 mt-1 rounded-lg bg-back text-left border border-back1 px-4"
            >
              Open <span className="hidden lg:inline-block">File</span>
              <input 
                onChange={(e) => {
                  const f = e.target.files[0];
                  setViewerLoaded(false);
                  const reader = new FileReader();
                  reader.readAsText(f);
                  setLoadingMessage("Loading " + f?.name)

                  reader.addEventListener(
                    "load",
                    () =>  { viewFile(reader.result, f.name); setOpenSelection(f.name)},
                    false
                  );
                }} 
                id="chemical-file" name="chemical-file" 
                type="file" accept=".pdb,.mol" 
                className="hidden" ref={fileUploadRef}
              />
            </PrimaryButton>
          </div>
        </div>

        <div className="relative flex flex-grow">
          <Modal
            isOpen={isHelpOpen}
            setIsOpen={setIsHelpOpen}
            className="max-w-md"
          >
            <h1 className="text-fore1 text-center text-2xl font-semibold">Help</h1>
            <div className="flex flex-col text-fore1 px-5">
              <div className="flex flex-shrink flex-row gap-8 items-center">
                <DragClickIcon className="w-7 fill-fore1"/>
                <span>Drag left mouse to orbit view</span>
              </div>

              <div className="flex flex-shrink flex-row gap-8 items-center">
                <DragClickIcon className="w-7 fill-fore1 -scale-x-100"/>
                <span>Drag right mouse to pan view</span>
              </div>

              <div className="flex flex-shrink flex-row gap-8 items-center">
                <ZoomIcon className="w-7 fill-fore1"/>
                <span>Scroll to zoom in and out view</span>
              </div>
            </div>
            <div className="flex justify-center mt-5">
              <PrimaryButton
                className="px-5 py-2 rounded-md"
                onClick={() => setIsHelpOpen(false)}
              >
                Close
              </PrimaryButton>
            </div>
          </Modal>

          <Modal
            isOpen={isSettingsOpen}
            setIsOpen={setIsSettingsOpen}
            className="max-w-lg overflow-visible"
          >
            <h1 className="text-fore1 text-center text-2xl font-semibold pb-2">Settings</h1>
            {viewerLoaded && 
            <div className="flex flex-col text-fore1 fill-fore gap-5">
              {(() => {
                const options = ["Perspective", "Orthographic"]
                const modes = {
                    [viewer.CameraType.PERSPECTIVE.value]: "Perspective",
                    [viewer.CameraType.ORTHOGRAPHIC.value]: "Orthographic",
                    "Perspective": viewer.CameraType.PERSPECTIVE,
                    "Orthographic": viewer.CameraType.ORTHOGRAPHIC,
                };

                return <div className='flex flex-row items-center'>
                    <span className='pr-4 w-[33%] min-w-[33%]'>Camera Projection</span>
                    <Combobox 
                        optionClassName="z-10"
                        className="flex-grow"
                        selected={modes[cameraType.value]}
                        setSelected={(opt) => setCameraType(modes[opt])}
                        options={options}
                    />
                </div>
              })()}

              {(() => {
                const options = ["Normal", "Gooch", "Simple"]
                const modes = {
                    [viewer.DrawSettingsMode.NORMAL.value]: "Normal",
                    [viewer.DrawSettingsMode.GOOCH.value]: "Gooch",
                    [viewer.DrawSettingsMode.BLINN_PHONG.value]: "Simple",
                    "Normal": viewer.DrawSettingsMode.NORMAL,
                    "Gooch": viewer.DrawSettingsMode.GOOCH,
                    "Simple": viewer.DrawSettingsMode.BLINN_PHONG,
                };

                return <div className='flex flex-row items-center'>
                    <span className='pr-4 w-[33%] min-w-[33%]'>Shading</span>
                    <Combobox 
                        optionClassName="z-10"
                        className="flex-grow"
                        selected={modes[drawSettings?.mode.value]}
                        setSelected={(opt) => setDrawSettings({...drawSettings, mode: modes[opt]})}
                        options={options}
                    />
                  </div>
                })()}
            </div>}
            <div className="flex justify-center mt-5">
              <PrimaryButton
                className="px-5 py-2 rounded-md"
                onClick={() => setIsSettingsOpen(false)}
              >
                Close
              </PrimaryButton>
            </div>
          </Modal>

          {/* Viewer canvas, must have id "canvas" */}
          {!viewerLoaded && 
          <>
            <div className="absolute bg-back w-full h-full z-10" />
            <div className="m-auto z-20 flex flex-col gap-6 items-center">
              <Throbber svgClassName="h-40 w-40"/>
              <span className="font-bold text-accent1 dark:text-accent">{loadingMessage}</span>
            </div> 
          </>
          }
          <canvas 
            onContextMenu={(e) => e.preventDefault()}
            className="absolute w-full h-full cursor-pointer" 
            id="canvas" ref={canvasRef}
          >
            This web browser doesn&lsquo;t support canvas tag.
          </canvas>
          {viewerLoaded && 
          <div className="absolute right-2 top-2 z-10 flex flex-col gap-3 items-center">
            <Tooltip padding={5} placement={"left"}>
              <TooltipTrigger>
                <PrimaryButton 
                  onClick={() => setIsSettingsOpen(!isSettingsOpen)}
                  className="rounded-full"
                >
                  <SettingsIcon className="w-10 h-10 p-2"/>
                </PrimaryButton>
              </TooltipTrigger>
              
              <TooltipContent>
                <div className="bg-back p-2 border border-back1 rounded-md">
                  <span className="block m-auto text-fore1 text-sm">
                    Open settings
                  </span>
                </div>
              </TooltipContent>
            </Tooltip>

            <Tooltip padding={5} placement={"left"}>
              <TooltipTrigger>
                <PrimaryButton 
                  onClick={() => setIsHelpOpen(!isHelpOpen)}
                  className="rounded-full"
                >
                  <HelpIcon className="w-10 h-10 p-2"/>
                </PrimaryButton>
              </TooltipTrigger>
              
              <TooltipContent>
                <div className="bg-back p-2 border border-back1 rounded-md">
                  <span className="block m-auto text-fore1 text-sm">
                    Open help
                  </span>
                </div>
              </TooltipContent>
            </Tooltip>

            <PrimaryButton 
              onClick={() => setFullRightPanel(!fullRightPanel)}
              className="lg:hidden rounded-full"
            >
              <MoreVertIcon className="w-10 h-10 p-2"/>
            </PrimaryButton>
          </div>
          }
        </div>
      </Panel>
      
      <PanelResizeHandle className="hidden lg:flex bg-back1 bg-opacity-10 w-3 hover:bg-opacity-20 hover:border-dashed">
        <button onClick={() => togglePanel(rightPanel)} className="max-w-full py-1 my-auto scale-[220%] group">
          <ChevronRightIcon className={
            classnames("w-[100%] transition-transform", {
              "rotate-180": rightCollapsed,
              "group-hover:translate-x-[2px]": !rightCollapsed,
              "group-hover:-translate-x-[2px]": rightCollapsed,
            })
          }/>
        </button>
      </PanelResizeHandle>

      <Panel className="hidden lg:flex border-l border-back1" defaultSize={30} minSize={20} collapsible={true} ref={rightPanel} onCollapse={setRightCollapsed}>
        <RightPanel viewer={viewer} isProtein={isProtein} refresh={refresh} viewerLoaded={viewerLoaded} setViewerLoaded={setViewerLoaded} setLoadingMessage={setLoadingMessage}/>
      </Panel>

      <Modal 
        isOpen={fullRightPanel}
        setIsOpen={setFullRightPanel}
        className="flex flex-col items-center bg-back z-10 border-back1 border-l p-0 text-fore fill-fore"
      >
        <RightPanel viewer={viewer} isProtein={isProtein} refresh={refresh} viewerLoaded={viewerLoaded} setViewerLoaded={setViewerLoaded} setLoadingMessage={setLoadingMessage}/>
        <PrimaryButton
          className="px-5 py-2 rounded-md mb-2 w-[10rem]"
          onClick={() => setFullRightPanel(false)}
        >
          Close
        </PrimaryButton>
      </Modal>
    </PanelGroup>
  )
}

export default App;
