import classnames from 'classnames';

import { ReactComponent as CheckIcon } from 'assets/check.svg';

const Check = ({ value, setValue, className }) => {
  return (
    <span
      className={classnames(
        'inline-block w-5 h-5 rounded-md hover:bg-transparent relative group hover:cursor-pointer bg-accent1 bg-opacity-40',
        {
          'bg-transparent': value,
          [className]: className
        }
      )}
      onClick={() => {
        setValue(!value);
      }}
    >
      <span className="w-8 h-8 -top-2 -left-1 absolute rounded-full group-hover:bg-back1"></span>

      <CheckIcon
        className={classnames('w-7 h-7 -top-1 absolute', {
          'invisible': !value,
          'visible': value
        })}
      ></CheckIcon>
    </span>
  );
};

export default Check;
