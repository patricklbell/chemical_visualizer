import { useEffect, useState } from "react";
import useExternalScript from "./useExternalScript";

export const useWasm = (canvasRef) => {
  const [module, setModule] = useState(null);
  const [script, setScript] = useState(undefined);
  useExternalScript(script);

  useEffect(
    (debug = true, print = true) => {
      if (!canvasRef || globalThis?.Module) return;

      const fetchWasm = async () => {
        fetch("./wasm/cviz_wasm.wasm")
          .then((response) => response.arrayBuffer())
          .then((bytes) => {
            globalThis.Module = {
              wasmBinary: bytes,
              preRun: [],
              postRun: [],
              canvas: canvasRef.current,
              // Once the wasm file has been loaded, the onRuntimeInitialized callback
              // will be invoked.
              onRuntimeInitialized: function () {
                setModule(globalThis.Module);
              },
              printErr: function () {
                if (debug)
                  console.warn(Array.prototype.slice.call(arguments).join(" "));
              },
              print: function () {
                if (print)
                  console.log(Array.prototype.slice.call(arguments).join(" "));
              },
            };

            setScript("./wasm/cviz_wasm.js");
          });
      };

      fetchWasm();
    },
    [canvasRef]
  );

  return module;
};

export const useSetting = (module, name) => {
  const [setting, setSetting] = useState(undefined);

  useEffect(() => {
    if (!setting) setSetting(module?.["get" + name]());
  }, [module]);

  useEffect(() => {
    if (!setting) return;

    module?.["set" + name](setting);
  }, [module, setting]);

  return [setting, setSetting];
};

// Converts all wasm vectors recursively to js list
const wasmToJs = (wsm) => {
  const vectorToJs = (vec) => {
    let arr = [];
    for (let i = 0; i < vec?.size(); i++) {
      arr.push(wasmToJs(vec.get(i)));
    }
    return arr;
  };
  const isVector = (obj) => {
    return !!(obj?.size && obj?.get);
  };

  if (!wsm || typeof wsm !== "object") return wsm;

  let obj = {};
  for (const prop in wsm) {
    const val = wsm[prop];

    obj[prop] = isVector(val) ? vectorToJs(val) : wasmToJs(val);
  }

  return obj !== {} ? obj : wsm;
};

export const useValue = (module, name, trigger) => {
  const [value, setValue] = useState(undefined);

  useEffect(() => {
    setValue(wasmToJs(module?.["get" + name]()));
  }, [module, trigger]);

  return [value];
};

export const usePdbDrawSettings = (module) =>
  useSetting(module, "PdbDrawSettings");

export const useDrawSettings = (module) => useSetting(module, "DrawSettings");

export const useCameraType = (module) => useSetting(module, "CameraType");

export const usePdbInfo = (module, refresh = true) =>
  useValue(module, "PdbInfo", refresh);

export const useMolInfo = (module, refresh = true) =>
  useValue(module, "MolInfo", refresh);
