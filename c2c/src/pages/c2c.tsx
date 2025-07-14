import Card from '../components/Card'

const C2C = () => {


  const copyToClipboard = (text: string) => {
    navigator.clipboard.writeText(text.replace(/\n/g, "\r\n"))
  };

  const downloadfile = (url: string) => {
    const fileName = url.split('/').pop();
    const a = document.createElement('a');
    a.href = url;
    a.setAttribute('download', fileName || 'download');
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }

  return (
    <div className='bg-[#0D1117] h-dvh w-auto flex justify-center items-center'>
      <div className='flex gap-5 pt-24 pb-10'>
        <Card classname="bg-[#584F44] text-[#9E9893]" text='Forza' onClick={() => downloadfile('./assets/Prob&StatDanveerDagur.pdf')} />
      </div>
    </div>
  );
}

export default C2C;
