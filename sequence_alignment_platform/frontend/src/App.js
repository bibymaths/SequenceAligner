import React, { useState, useEffect } from 'react';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import DropZone from './components/DropZone';
import AlignmentViewer from './components/AlignmentViewer';
import MatrixVisualizer from './components/MatrixVisualizer';
import AnalysisDashboard from './components/AnalysisDashboard';

const queryClient = new QueryClient();

function App() {
    // 1. State to manage the app's current phase
    const [sessionId, setSessionId] = useState(null);
    const [logs, setLogs] = useState([]);
    const [isComplete, setIsComplete] = useState(false);

    // 2. The WebSocket Listener
    useEffect(() => {
        // Don't connect if we don't have a session ID yet, or if we are already done
        if (!sessionId || isComplete) return;

        const ws = new WebSocket(`ws://127.0.0.1:8000/ws/logs/${sessionId}`);

        ws.onmessage = (event) => {
            const message = event.data;

            // Append new logs to the console
            setLogs((prevLogs) => [...prevLogs, message]);

            // THE MAGIC TRIGGER: If the C++ backend says it's done, flip the switch!
            if (message.includes("Analysis complete") || message.includes("completed successfully")) {
                setIsComplete(true);
                ws.close(); // Hang up the connection
            }
        };

        // Cleanup function if the component unmounts
        return () => {
            if (ws.readyState === 1) ws.close();
        };
    }, [sessionId, isComplete]);

    return (
        <QueryClientProvider client={queryClient}>
            <div className="min-h-screen bg-gray-100 p-8">
                <h1 className="text-3xl font-bold mb-6 text-gray-800">Sequence Alignment Platform</h1>

                {/* PHASE 1: Upload (Only show if no session exists) */}
                {!sessionId && (
                    <DropZone onSessionCreated={(id) => setSessionId(id)} />
                )}

                {/* PHASE 2: Running & Streaming Logs (Show if we have an ID but aren't done) */}
                {sessionId && !isComplete && (
                    <div className="bg-white p-6 rounded-lg shadow-md mt-4 border border-gray-200">
                        <h2 className="text-xl font-semibold mb-4 text-blue-600">
                            <span className="animate-pulse mr-2">⚙️</span>
                            Running C++ Alignment...
                        </h2>
                        <div className="bg-gray-900 text-green-400 p-4 rounded h-64 overflow-y-auto font-mono text-sm whitespace-pre-wrap">
                            {logs.map((log, index) => (
                                <span key={index}>{log}</span>
                            ))}
                        </div>
                    </div>
                )}

                {/* PHASE 3: Results Dashboard (Show ONLY when complete) */}
                {isComplete && (
                    <div className="space-y-8 mt-6">
                        <div className="bg-green-100 border-l-4 border-green-500 text-green-700 p-4 rounded shadow-sm">
                            <p className="font-bold">Success!</p>
                            <p>Alignment and analysis are complete. Rendering results below.</p>
                        </div>

                        {/* Pass the sessionId as a prop so they know which data to fetch! */}
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