import React, { useState } from 'react';
import axios from 'axios';

function DropZone({ onSessionCreated }) {
  // File states
  const [queryFile, setQueryFile] = useState(null);
  const [targetFile, setTargetFile] = useState(null);

  // Alignment parameter states
  const [seqType, setSeqType] = useState('dna');
  const [alignMethod, setAlignMethod] = useState('global');

  // UI states
  const [isUploading, setIsUploading] = useState(false);
  const [error, setError] = useState('');

  const handleUpload = async () => {
    if (!queryFile || !targetFile) {
      setError('Please select both a Query and Target file.');
      return;
    }

    setIsUploading(true);
    setError('');

    const formData = new FormData();
    formData.append('query', queryFile);
    formData.append('target', targetFile);
    formData.append('seq_type', seqType);           // 'dna' or 'protein'
    formData.append('align_method', alignMethod);   // 'global', 'local', etc.

    try {
      const response = await axios.post('http://127.0.0.1:8000/align', formData);

      // Pass the ID up to App.js to switch the view to the Live Logs
      if (onSessionCreated) {
        onSessionCreated(response.data.session_id);
      }
    } catch (err) {
      console.error("Upload error:", err);
      setError(
          err.response?.data?.detail || 'Failed to connect to the backend. Is FastAPI running?'
      );
    } finally {
      setIsUploading(false);
    }
  };

  return (
      <div className="bg-white p-6 rounded-lg shadow-md border border-gray-200 max-w-3xl">
        <h2 className="text-xl font-bold mb-4 text-gray-800">New Alignment Session</h2>

        {/* File Inputs */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
          <div className="border-2 border-dashed border-gray-300 rounded-md p-4 text-center hover:bg-gray-50 transition-colors">
            <label className="block text-sm font-medium text-gray-700 mb-2">Query Sequence (.fasta)</label>
            <input
                type="file"
                accept=".fa,.fasta"
                onChange={(e) => setQueryFile(e.target.files[0])}
                className="text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-md file:border-0 file:text-sm file:font-semibold file:bg-blue-50 file:text-blue-700 hover:file:bg-blue-100 cursor-pointer"
            />
          </div>

          <div className="border-2 border-dashed border-gray-300 rounded-md p-4 text-center hover:bg-gray-50 transition-colors">
            <label className="block text-sm font-medium text-gray-700 mb-2">Target Sequence (.fasta)</label>
            <input
                type="file"
                accept=".fa,.fasta"
                onChange={(e) => setTargetFile(e.target.files[0])}
                className="text-sm text-gray-500 file:mr-4 file:py-2 file:px-4 file:rounded-md file:border-0 file:text-sm file:font-semibold file:bg-blue-50 file:text-blue-700 hover:file:bg-blue-100 cursor-pointer"
            />
          </div>
        </div>

        {/* Alignment Parameters */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
          <div>
            <label className="block text-xs font-semibold text-gray-600 uppercase tracking-wide mb-1">Type</label>
            <select
                value={seqType}
                onChange={(e) => setSeqType(e.target.value)}
                className="w-full bg-gray-50 border border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2.5"
            >
              <option value="dna">DNA</option>
              <option value="protein">Protein</option>
            </select>
          </div>

          <div>
            <label className="block text-xs font-semibold text-gray-600 uppercase tracking-wide mb-1">Algorithm</label>
            <select
                value={alignMethod}
                onChange={(e) => setAlignMethod(e.target.value)}
                className="w-full bg-gray-50 border border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2.5"
            >
              <option value="global">Global (Needleman-Wunsch)</option>
              <option value="local">Local (Smith-Waterman)</option>
              <option value="lcs">Longest Common Subsequence</option>
              <option value="all">Run All Three</option>
            </select>
          </div>
        </div>

        {/* Error Message */}
        {error && (
            <div className="mb-4 text-sm text-red-600 bg-red-50 border border-red-200 p-3 rounded-md">
              {error}
            </div>
        )}

        {/* Submit Button */}
        <button
            onClick={handleUpload}
            disabled={isUploading}
            className={`w-full font-bold py-3 px-4 rounded-lg transition-colors ${
                isUploading
                    ? 'bg-blue-300 cursor-not-allowed text-white'
                    : 'bg-blue-600 hover:bg-blue-700 text-white shadow-md'
            }`}
        >
          {isUploading ? 'Initializing Alignment...' : 'Run Alignment'}
        </button>
      </div>
  );
}

export default DropZone;