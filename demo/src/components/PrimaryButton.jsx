import classnames from 'classnames';
import { useId } from 'react';

// Based on https://tailwind-elements.com/docs/standard/components/buttons/
const PrimaryButton = ({
  type = 'button',
  children,
  className,
  ...other
}) => {
  const id = useId();

  return (
    <button
      type={type}
      id={id}
      className={classnames(
        'text-sm border border-back1 bg-back text-fore1 transition duration-150 ease-in-out hover:bg-accent3 hover:border-accent1 focus:outline-none focus:ring-0 disabled:bg-back1 disabled:border-back',
        {
          [className]: className
        }
      )}
      {...other}
    >
      {children}
    </button>
  );
};

export default PrimaryButton;
