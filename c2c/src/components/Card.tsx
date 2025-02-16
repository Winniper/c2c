import clsx from "clsx"

type CardProps = {
    classname?: string;
    text: string;
    onClick?: () => void;
  };
const cardStyles = "flex items-center font-poppins font-bold text-sm sm:text-2xl md:text-3xl transition-all shadow-[3px_3px_0px_black] hover:shadow-none hover:translate-x-[5px] hover:translate-y-[5px] active:scale-95 p-0 sm:p-6 md:p-8 rounded-3xl"

function Card({classname, text, onClick}: CardProps) {
  return (
    <div onClick={onClick} className={clsx(cardStyles, classname)}>
        <span className="mx-2">{text}</span>
    </div>
  )
}

export default Card