import React, { useState, useEffect } from 'react';
import axios from 'axios';

function AnalysisDashboard({ sessionId }) {
    const [analysisData, setAnalysisData] = useState(null);
    const [tableData, setTableData] = useState({});
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState('');

    useEffect(() => {
        if (!sessionId) return;

        const fetchDashboardData = async () => {
            try {
                setLoading(true);
                // 1. Discover what analysis files exist
                const analysisRes = await axios.get(`http://127.0.0.1:8000/session/${sessionId}/analysis`);
                const groupedFiles = analysisRes.data;
                setAnalysisData(groupedFiles);

                // 2. Fetch the actual content of the TSV tables
                const fetchedTables = {};
                for (const group of Object.keys(groupedFiles)) {
                    const tsvFiles = groupedFiles[group].tsv || [];
                    for (const filename of tsvFiles) {
                        const tableRes = await axios.get(`http://127.0.0.1:8000/session/${sessionId}/analysis/table/${filename}`);
                        fetchedTables[filename] = tableRes.data.records;
                    }
                }
                setTableData(fetchedTables);

            } catch (err) {
                console.error("Failed to load analysis dashboard:", err);
                setError("Analysis data not found. Downstream analysis may have failed or is still running.");
            } finally {
                setLoading(false);
            }
        };

        fetchDashboardData();
    }, [sessionId]);

    if (!sessionId) return null;

    return (
        <div className="bg-white rounded-lg shadow-md border border-gray-200 mt-6 p-6">
            <h3 className="text-xl font-bold mb-4 text-gray-800 border-b pb-2">Analysis Dashboard</h3>

            {loading && <p className="text-blue-500 animate-pulse font-semibold">Loading analytical data...</p>}
            {error && <p className="text-red-500 bg-red-50 p-3 rounded">{error}</p>}

            {!loading && !error && analysisData && (
                <div className="space-y-8 mt-4">
                    {Object.keys(analysisData).map((groupKey) => (
                        <div key={groupKey} className="bg-gray-50 border rounded-lg p-4 shadow-sm">
                            <h4 className="text-lg font-bold text-blue-800 uppercase tracking-wider mb-4">
                                {groupKey} Analysis
                            </h4>

                            {/* Render TSV Tables */}
                            {analysisData[groupKey].tsv?.map((filename) => (
                                <div key={filename} className="mb-6 overflow-x-auto">
                                    <h5 className="text-sm font-semibold text-gray-600 mb-2">{filename}</h5>
                                    <table className="min-w-full bg-white border border-gray-200 rounded text-sm text-left">
                                        <thead className="bg-gray-100 text-gray-700">
                                        <tr>
                                            {tableData[filename] && tableData[filename].length > 0 &&
                                                Object.keys(tableData[filename][0]).map((header) => (
                                                    <th key={header} className="py-2 px-4 border-b font-semibold">{header}</th>
                                                ))
                                            }
                                        </tr>
                                        </thead>
                                        <tbody>
                                        {tableData[filename]?.map((row, idx) => (
                                            <tr key={idx} className="hover:bg-blue-50 transition-colors">
                                                {Object.values(row).map((val, vIdx) => (
                                                    <td key={vIdx} className="py-2 px-4 border-b text-gray-700">{val}</td>
                                                ))}
                                            </tr>
                                        ))}
                                        </tbody>
                                    </table>
                                </div>
                            ))}

                            {/* Render PNG Images */}
                            {analysisData[groupKey].png?.length > 0 && (
                                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                                    {analysisData[groupKey].png.map((filename) => (
                                        <div key={filename} className="bg-white border p-2 rounded shadow-sm text-center">
                                            <h5 className="text-xs font-semibold text-gray-500 mb-2">{filename}</h5>
                                            <img
                                                src={`http://127.0.0.1:8000/session/${sessionId}/file/analysis_out/${filename}`}
                                                alt={filename}
                                                className="max-w-full h-auto mx-auto rounded"
                                            />
                                        </div>
                                    ))}
                                </div>
                            )}
                        </div>
                    ))}
                </div>
            )}
        </div>
    );
}

export default AnalysisDashboard;