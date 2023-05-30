import { Transition, Dialog } from '@headlessui/react';
import { Fragment } from 'react';
import classnames from 'classnames';

const Modal = (props) => {
  const { isOpen, setIsOpen } = props;
  const { closeModal = () => setIsOpen(false) } = props;

  // Modified from https://headlessui.com/react/dialog
  return (
    <Transition appear show={isOpen} as={Fragment}>
      <Dialog as="div" className="relative z-10" onClose={closeModal}>
        <Transition.Child
          as={Fragment}
          enter="ease-out duration-200"
          enterFrom="opacity-0"
          enterTo="opacity-100"
          leave="ease-in duration-200"
          leaveFrom="opacity-100"
          leaveTo="opacity-0"
        >
          <div className="fixed inset-0 bg-fore1 bg-opacity-25" />
        </Transition.Child>

        <div className="fixed inset-0 overflow-y-auto">
          <div className="flex min-h-full items-center justify-center p-4 text-center">
            <Transition.Child
              as={Fragment}
              enter="ease-out duration-200"
              enterFrom="opacity-0 scale-95"
              enterTo="opacity-100 scale-100"
              leave="ease-in duration-200"
              leaveFrom="opacity-100 scale-100"
              leaveTo="opacity-0 scale-95"
            >
              <Dialog.Panel
                className={classnames(
                  'w-full max-w-2xl transform overflow-hidden rounded-2xl bg-back p-6 text-left align-middle shadow-xl transition-all',
                  {
                    [props.className]: props.className
                  }
                )}
              >
                <Dialog.Title
                  as="h3"
                  className="text-lg font-medium leading-6 text-fore1"
                >
                  {props?.title}
                </Dialog.Title>
                {props?.children}
              </Dialog.Panel>
            </Transition.Child>
          </div>
        </div>
      </Dialog>
    </Transition>
  );
};

export default Modal;
