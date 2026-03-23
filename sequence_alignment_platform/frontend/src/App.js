import React, { useState, useEffect } from 'react';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import DropZone from './components/DropZone';
import AlignmentViewer from './components/AlignmentViewer';
import MatrixVisualizer from './components/MatrixVisualizer';
import AnalysisDashboard from './components/AnalysisDashboard';

const queryClient = new QueryClient();

function App() {
    const [sessionId, setSessionId] = useState(null);
    const [status, setStatus] = useState("idle"); // idle | running | completed | failed

    useEffect(() => {
        if (!sessionId || status !== "running") return;

        const ws = new WebSocket(`ws://127.0.0.1:8000/ws/logs/${sessionId}`);

        ws.onmessage = (event) => {
            const message = event.data.toLowerCase();

            if (message.includes("[error]")) {
                setStatus("failed");
                ws.close();
                return;
            }

            if (
                message.includes("analysis complete") ||
                message.includes("alignment complete") ||
                message.includes("pipeline finished successfully")
            ) {
                setStatus("completed");
                ws.close();
            }
        };

        ws.onerror = () => {
            setStatus("failed");
        };

        return () => {
            if (ws.readyState === WebSocket.OPEN || ws.readyState === WebSocket.CONNECTING) {
                ws.close();
            }
        };
    }, [sessionId, status]);

    const resetRun = () => {
        setSessionId(null);
        setStatus("idle");
    };

    return (
        <QueryClientProvider client={queryClient}>
            <div className="min-h-screen bg-gray-100 p-8">
                <h1 className="text-3xl font-bold mb-6 text-gray-800">
                    Sequence Alignment Platform
                </h1>

                {/* Upload */}
                {!sessionId && (
                    <DropZone
                        onSessionCreated={(id) => {
                            setSessionId(id);
                            setStatus("running");
                        }}
                    />
                )}

                {/* Running State */}
                {sessionId && status === "running" && (
                    <div className="bg-white p-6 rounded-2xl shadow-md border border-gray-200 mt-4">
                        <div className="flex items-center gap-4">
                            <div className="animate-spin h-5 w-5 border-2 border-blue-500 border-t-transparent rounded-full"></div>
                            <div>
                                <h2 className="text-lg font-semibold text-gray-800">
                                    Running alignment & analysis
                                </h2>
                                <p className="text-sm text-gray-500">
                                    This may take a few seconds depending on sequence size.
                                </p>
                            </div>
                        </div>
                    </div>
                )}

                {/* Failed */}
                {sessionId && status === "failed" && (
                    <div className="bg-red-100 border-l-4 border-red-500 text-red-700 p-4 rounded shadow-sm mt-6">
                        <p className="font-bold">Pipeline failed</p>
                        <p>Check backend logs. This UI intentionally hides raw terminal output.</p>

                        <button
                            onClick={resetRun}
                            className="mt-4 px-3 py-2 rounded bg-gray-200 hover:bg-gray-300 text-sm"
                        >
                            New Run
                        </button>
                    </div>
                )}

                {/* Success */}
                {sessionId && status === "completed" && (
                    <div className="space-y-8 mt-6">
                        <div className="bg-green-100 border-l-4 border-green-500 text-green-700 p-4 rounded shadow-sm flex justify-between items-center">
                            <div>
                                <p className="font-bold">Completed</p>
                                <p>Alignment and analysis finished successfully.</p>
                            </div>

                            <button
                                onClick={resetRun}
                                className="px-3 py-2 rounded bg-gray-200 hover:bg-gray-300 text-sm"
                            >
                                New Run
                            </button>
                        </div>

                        <AlignmentViewer sessionId={sessionId} />
                        <MatrixVisualizer sessionId={sessionId} />
                        <AnalysisDashboard sessionId={sessionId} />
                    </div>
                )}
            </div>
        </QueryClientProvider>
    );
}

export default App;