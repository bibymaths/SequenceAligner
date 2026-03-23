import React, { useState } from 'react';
import { FixedSizeList as List } from 'react-window';

// Placeholder component for streaming large alignments. In a real implementation,
// this would fetch segments of the alignment from the backend and render
// colored spans for matches, mismatches, and gaps. It uses react-window
// to virtualize rows of alignment, ensuring that only visible lines are
// rendered to the DOM at any time.

const dummyData = Array.from({ length: 10000 }, (_, i) => {
  return {
    index: i,
    seq1: 'ACGTACGT'.repeat(10),
    seq2: 'ACGTTCGT'.repeat(10),
  };
});

const Row = ({ index, style }) => {
  const row = dummyData[index];
  return (
    <div style={style} className="font-mono whitespace-pre">
      {row.seq1}
      <br />
      {row.seq2}
    </div>
  );
};

function AlignmentViewer() {
  return (
    <div className="border p-2 mt-4 bg-white" style={{ height: '300px' }}>
      <h3 className="font-semibold mb-2">Alignment Viewer</h3>
      <List height={250} itemCount={dummyData.length} itemSize={32} width="100%">
        {Row}
      </List>
    </div>
  );
}

export default AlignmentViewer;