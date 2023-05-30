import { Disclosure, Transition } from '@headlessui/react';
import classnames from 'classnames';

import { ReactComponent as ExpandIcon } from 'assets/Expand.svg';

export const Accordian = ({ children, ...props }) => 
    <Disclosure {...props}>
      {children}
    </Disclosure>

export const AccordianButton = ({ children, className, open, ...props }) => 
    <Disclosure.Button 
      className={classnames('bg-back2 w-full text-left p-1 m-1 ml-0 rounded-md hover:bg-accent3', {[className]: className})}
      {...props}
    >
    {({ open }) => (
      <div className='flex flex-row gap-1 items-center'>
        <ExpandIcon className={classnames('w-7 h-7 pt-1 transition-transform', {'-scale-y-100': open})}/>
        <div>{children}</div>
      </div>
    )}
    </Disclosure.Button>

export const AccordianPanel = ({ children, className, ...props }) => 
    <Transition
        enter="transition duration-100 ease-out"
        enterFrom="transform scale-95 opacity-0"
        enterTo="transform scale-100 opacity-100"
        leave="transition duration-75 ease-out"
        leaveFrom="transform scale-100 opacity-100"
        leaveTo="transform scale-95 opacity-0"
    >

      <Disclosure.Panel 
        className={classnames('pl-8 py-1 my-1', {[className]: className})}
        {...props}
      >
        {children}
      </Disclosure.Panel>
    </Transition>
