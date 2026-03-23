import React, { useState } from 'react';
import axios from 'axios';

function DropZone() {
  const [queryFile, setQueryFile] = useState(null);
  const [targetFile, setTargetFile] = useState(null);
  const [status, setStatus] = useState('');

  const handleUpload = async () => {
    if (!queryFile || !targetFile) return;
    const formData = new FormData();
    formData.append('query', queryFile);
    formData.append('target', targetFile);
    formData.append('gap_open', 10);
    formData.append('gap_extend', 0.5);
    formData.append('mode', 'global');
    try {
      const res = await axios.post('/align', formData);
      setStatus(`Session created: ${res.data.session_id}`);
    } catch (err) {
      console.error(err);
      setStatus('Upload failed');
    }
  };

  return (
    <div className="border-dashed border-2 p-4 mb-4 bg-white">
      <h2 className="font-semibold mb-2">Drop & Align</h2>
      <div className="flex space-x-4 mb-2">
        <input type="file" accept=".fa,.fasta" onChange={(e) => setQueryFile(e.target.files[0])} />
        <input type="file" accept=".fa,.fasta" onChange={(e) => setTargetFile(e.target.files[0])} />
      </div>
      <button className="bg-blue-500 text-white px-4 py-2" onClick={handleUpload}>Run Alignment</button>
      {status && <p className="mt-2 text-sm text-gray-600">{status}</p>}
    </div>
  );
}

export default DropZone;