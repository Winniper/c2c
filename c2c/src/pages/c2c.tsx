import Card from '../components/Card'

const C2C = () => {
  const onClickGraph = () =>{
    console.log('graph')
  }

  const onClickBFS = () => {
    console.log('BFS')
  }

  return (
    <div className='bg-[#0D1117] h-dvh w-auto flex justify-center items-center'>
      <div className='flex flex-col gap-5 w-11/12 sm:w-3/4 h-full pt-24 pb-10'>
        <Card classname="bg-[#29903B] text-white" text='Graphical' onClick={onClickGraph} />
        <Card classname="bg-[#29903B] text-white" text='BFS' onClick={onClickBFS}/>
      </div>
    </div>
  )
}

export default C2C