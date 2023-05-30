import { Combobox, Transition } from '@headlessui/react';
import { Fragment, useEffect, useState } from 'react';
import classnames from 'classnames';

import { ReactComponent as ChevronRightIcon } from 'assets/ChevronRight.svg';
import { ReactComponent as CheckIcon } from 'assets/Check.svg';

const InputCombobox = ({
  selected = '',
  className = '',
  optionClassName = '',
  inputClassName = '',
  setSelected,
  options,
  placeholder,
  filter = (option, query) =>
    option
      .toLowerCase()
      .replace(/\s+/g, '')
      .includes(query.toLowerCase().replace(/\s+/g, '')),
  setInput = () => {}
}) => {
  const [query, setQuery] = useState('');
  useEffect(() => setInput(query), [query]);

  const filteredOptions = options?.filter
    ? options.filter((option) => filter(option, query))
    : options;

  return (
    <Combobox value={selected} onChange={setSelected} className={classnames("group", {[className]: className})}>
      <div className="relative mt-1">
        <div className="relative w-full cursor-default overflow-hidden rounded-lg bg-back text-left border border-back1">
          <Combobox.Input
            className={classnames(
              "w-full border-none py-3 pl-3 pr-10 text-sm leading-5 bg-back text-fore1", 
              {[inputClassName]: inputClassName}
            )}
            displayValue={(option) => option}
            onChange={(event) => setQuery(event.target.value)}
            placeholder={placeholder}
          />
          <Combobox.Button
            className={classnames(
              'absolute inset-y-0 right-0 flex items-center pr-2',
              { invisible: !options?.length }
            )}
          >
            <ChevronRightIcon
              className="h-5 w-5 text-fore1 rotate-90"
              aria-hidden="true"
            />
          </Combobox.Button>
        </div>
        <Transition
          as={Fragment}
          className="z-10"
          leave="transition ease-in duration-100"
          leaveFrom="opacity-100"
          leaveTo="opacity-0"
          afterLeave={() => setQuery('')}
        >
          <Combobox.Options
            className={classnames(
              'absolute mt-1 max-h-60 w-full overflow-auto rounded-md bg-back text-fore border border-back1 shadow-md py-1 text-base focus:outline-none sm:text-sm',
              { invisible: !options?.length }
            )}
          >
            {filteredOptions.length === 0 && query !== '' ? (
              <div className="relative cursor-default select-none py-2 px-4">
                Nothing found.
              </div>
            ) : (
              filteredOptions.map((option) => (
                <Combobox.Option
                  tabindex="-1"
                  key={option}
                  className={({ active }) =>
                    classnames(
                      'relative cursor-pointer select-none py-2 pl-10 pr-4',
                      {'bg-accent2 bg-opacity-20': active, [optionClassName]: optionClassName }
                    )
                  }
                  value={option}
                >
                  {({ selected, active }) => (
                    <>
                      <span
                        className={classnames(
                          'block truncate', 
                          { 'text-accent font-bold' : selected }
                        )}
                      >
                        {option}
                      </span>
                      {selected ? (
                        <span
                          className='absolute inset-y-0 left-0 flex items-center pl-3'
                        >
                          <CheckIcon
                            className="h-5 w-5 fill-accent"
                            aria-hidden="true"
                          />
                        </span>
                      ) : null}
                    </>
                  )}
                </Combobox.Option>
              ))
            )}
          </Combobox.Options>
        </Transition>
      </div>
    </Combobox>
  );
};

export default InputCombobox;
