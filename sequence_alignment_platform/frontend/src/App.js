import React from 'react';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import DropZone from './components/DropZone';
import AlignmentViewer from './components/AlignmentViewer';
import MatrixVisualizer from './components/MatrixVisualizer';
import AnalysisDashboard from './components/AnalysisDashboard';

const queryClient = new QueryClient();

function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <div className="min-h-screen bg-gray-100 p-4">
        <h1 className="text-2xl font-bold mb-4">Sequence Alignment Platform</h1>
        <DropZone />
        {/* Additional components like live console, matrix visualizer, and analysis dashboard would be conditionally rendered based on session state */}
        <AlignmentViewer />
        <MatrixVisualizer />
        <AnalysisDashboard />
      </div>
    </QueryClientProvider>
  );
}

export default App;