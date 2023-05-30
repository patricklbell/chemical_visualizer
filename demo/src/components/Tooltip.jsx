// Based on https://floating-ui.com/docs/tooltip
import {
  useFloating,
  autoUpdate,
  offset,
  flip,
  shift,
  useHover,
  useFocus,
  useDismiss,
  useRole,
  useInteractions,
  useMergeRefs,
  FloatingPortal
} from '@floating-ui/react';
import {
  cloneElement,
  isValidElement,
  createContext,
  forwardRef,
  useContext,
  useMemo,
  useState
} from 'react';
import { Transition } from '@headlessui/react';

export function useTooltip({
  initialOpen = false,
  placement = 'bottom',
  padding = 0,
  doHover = true,
  open: controlledOpen,
  onOpenChange: setControlledOpen
}) {
  const [uncontrolledOpen, setUncontrolledOpen] = useState(initialOpen);

  const open = controlledOpen ?? uncontrolledOpen;
  const setOpen = setControlledOpen ?? setUncontrolledOpen;

  const data = useFloating({
    placement,
    open,
    onOpenChange: setOpen,
    whileElementsMounted: autoUpdate,
    middleware: [
      offset(padding),
      flip({
        fallbackAxisSideDirection: 'start'
      }),
      shift({ padding })
    ]
  });

  const context = data.context;

  const hover = useHover(context, {
    move: false,
    enabled: doHover && controlledOpen == null
  });
  const focus = useFocus(context, {
    enabled: controlledOpen == null
  });
  const dismiss = useDismiss(context);
  const role = useRole(context, { role: 'tooltip' });

  const interactions = useInteractions([hover, focus, dismiss, role]);

  return useMemo(
    () => ({
      open,
      setOpen,
      ...interactions,
      ...data
    }),
    [open, setOpen, interactions, data]
  );
}

const TooltipContext = createContext(null);

const useTooltipContext = () => {
  const context = useContext(TooltipContext);

  if (context == null) {
    throw new Error('Tooltip components must be wrapped in <Tooltip />');
  }

  return context;
};

export function Tooltip({ children, ...options }) {
  // This can accept any props as options, e.g. `placement`,
  // or other positioning options.
  const tooltip = useTooltip(options);
  return (
    <TooltipContext.Provider value={tooltip}>
      {children}
    </TooltipContext.Provider>
  );
}

export const TooltipTrigger = forwardRef(function TooltipTrigger(
  { children, asChild = false, ...props },
  propRef
) {
  const context = useTooltipContext();
  const childrenRef = children.ref;
  const ref = useMergeRefs([context.refs.setReference, propRef, childrenRef]);

  // `asChild` allows the user to pass any element as the anchor
  if (asChild && isValidElement(children)) {
    return cloneElement(
      children,
      context.getReferenceProps({
        ref,
        ...props,
        ...children.props,
        'data-state': context.open ? 'open' : 'closed'
      })
    );
  }

  return (
    <button
      type="button"
      tabIndex="-1"
      ref={ref}
      // The user can style the trigger based on the state
      data-state={context.open ? 'open' : 'closed'}
      {...context.getReferenceProps(props)}
    >
      {children}
    </button>
  );
});

export const TooltipContent = forwardRef(function TooltipContent(
  props,
  propRef
) {
  const context = useTooltipContext();
  const ref = useMergeRefs([context.refs.setFloating, propRef]);

  return (
    <FloatingPortal>
      <Transition
        show={context.open}
        appear={true}
        enter="transition duration-100 ease-out"
        enterFrom="opacity-0"
        enterTo="opacity-100"
        leave="transition duration-75 ease-out"
        leaveFrom="opacity-100"
        leaveTo="opacity-0"
      >
        <div
          ref={ref}
          style={{
            position: context.strategy,
            top: context.y ?? 0,
            left: context.x ?? 0,
            ...props.style
          }}
          {...context.getFloatingProps(props)}
        />
      </Transition>
    </FloatingPortal>
  );
});
