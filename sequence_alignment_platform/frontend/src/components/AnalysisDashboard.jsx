import React, { useState, useEffect, useMemo } from 'react';
import axios from 'axios';

const METHOD_ORDER = ['global', 'local', 'lcs', 'comparison', 'summary'];
const SECTION_ORDER = [
    'summary',
    'alignment_summary',
    'path_metrics',
    'conserved_blocks',
    'substitution_summary',
    'residue_support',
    'method_comparison',
    'method_comparison_categories',
    'heatmaps',
    'other',
];

function prettifyLabel(raw) {
    return raw
        .replace(/_/g, ' ')
        .replace(/\btsv\b|\bpng\b|\bjson\b/gi, '')
        .replace(/\s+/g, ' ')
        .trim()
        .replace(/\b\w/g, (c) => c.toUpperCase());
}

function detectMethod(filename) {
    const lower = filename.toLowerCase();
    if (lower.includes('_global_')) return 'global';
    if (lower.includes('_local_')) return 'local';
    if (lower.includes('_lcs_')) return 'lcs';
    if (lower.includes('alignment_method_comparison')) return 'comparison';
    if (lower.endsWith('_summary.json')) return 'summary';
    return 'other';
}

function detectSection(filename) {
    const lower = filename.toLowerCase();

    if (lower.endsWith('_summary.json')) return 'summary';
    if (lower.includes('alignment_summary')) return 'alignment_summary';
    if (lower.includes('path_metrics')) return 'path_metrics';
    if (lower.includes('conserved_blocks')) return 'conserved_blocks';
    if (lower.includes('substitution_summary')) return 'substitution_summary';
    if (lower.includes('residue_support')) return 'residue_support';
    if (lower.includes('alignment_method_comparison_categories')) return 'method_comparison_categories';
    if (lower.includes('alignment_method_comparison')) return 'method_comparison';
    if (lower.includes('dp_heatmap')) return 'heatmaps';

    return 'other';
}

function cleanTitle(filename, sessionId) {
    let name = filename;

    if (sessionId) {
        name = name.replace(`${sessionId}_`, '');
    }

    name = name.replace(/\.(tsv|png|json)$/i, '');

    name = name
        .replace(/_with_path/gi, ' with path')
        .replace(/_aligned/gi, ' aligned')
        .replace(/_/g, ' ')
        .replace(/\s+/g, ' ')
        .trim();

    return name.replace(/\b\w/g, (c) => c.toUpperCase());
}

function sortByOrder(items, extractor, order) {
    const indexMap = new Map(order.map((key, idx) => [key, idx]));
    return [...items].sort((a, b) => {
        const aKey = extractor(a);
        const bKey = extractor(b);
        const aIdx = indexMap.has(aKey) ? indexMap.get(aKey) : 999;
        const bIdx = indexMap.has(bKey) ? indexMap.get(bKey) : 999;
        if (aIdx !== bIdx) return aIdx - bIdx;
        return String(a).localeCompare(String(b));
    });
}

function SectionCard({ title, defaultOpen = true, children, rightNote = '' }) {
    const [open, setOpen] = useState(defaultOpen);

    return (
        <div className="border border-gray-200 rounded-2xl bg-white shadow-sm">
            <button
                type="button"
                onClick={() => setOpen(!open)}
                className="w-full flex items-center justify-between px-5 py-4 text-left"
            >
                <div>
                    <h4 className="text-base font-semibold text-gray-800">{title}</h4>
                    {rightNote ? (
                        <p className="text-xs text-gray-500 mt-1">{rightNote}</p>
                    ) : null}
                </div>
                <span className="text-sm text-gray-500">{open ? 'Hide' : 'Show'}</span>
            </button>

            {open && <div className="px-5 pb-5">{children}</div>}
        </div>
    );
}

function DataTable({ title, rows }) {
    if (!rows || rows.length === 0) {
        return (
            <div className="border border-dashed border-gray-300 rounded-xl p-4 text-sm text-gray-500 bg-gray-50">
                No rows available for {title}.
            </div>
        );
    }

    const headers = Object.keys(rows[0]);

    return (
        <div className="border border-gray-200 rounded-xl overflow-hidden bg-white">
            <div className="px-4 py-3 border-b bg-gray-50">
                <h5 className="text-sm font-semibold text-gray-700">{title}</h5>
                <p className="text-xs text-gray-500 mt-1">{rows.length} rows</p>
            </div>

            <div className="overflow-auto max-h-[420px]">
                <table className="min-w-full text-sm">
                    <thead className="bg-gray-100 sticky top-0 z-10">
                    <tr>
                        {headers.map((header) => (
                            <th
                                key={header}
                                className="px-4 py-3 text-left font-semibold text-gray-700 border-b whitespace-nowrap"
                            >
                                {prettifyLabel(header)}
                            </th>
                        ))}
                    </tr>
                    </thead>
                    <tbody>
                    {rows.map((row, idx) => (
                        <tr key={idx} className="odd:bg-white even:bg-gray-50 hover:bg-blue-50">
                            {headers.map((header) => (
                                <td
                                    key={`${idx}-${header}`}
                                    className="px-4 py-2 border-b text-gray-700 align-top whitespace-nowrap"
                                >
                                    {String(row[header] ?? '')}
                                </td>
                            ))}
                        </tr>
                    ))}
                    </tbody>
                </table>
            </div>
        </div>
    );
}

function ImageCard({ sessionId, filename, sessionLabel }) {
    return (
        <div className="border border-gray-200 rounded-2xl bg-white shadow-sm overflow-hidden">
            <div className="px-4 py-3 border-b bg-gray-50">
                <h5 className="text-sm font-semibold text-gray-700">
                    {cleanTitle(filename, sessionLabel)}
                </h5>
            </div>

            <div className="p-4 flex items-center justify-center bg-gray-50">
                <img
                    src={`/session/${sessionId}/file/analysis_out/${filename}`}
                    alt={filename}
                    className="w-full max-h-[420px] object-contain rounded-lg border border-gray-200 bg-white"
                    loading="lazy"
                />
            </div>
        </div>
    );
}

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
                setError('');

                const analysisRes = await axios.get(`/session/${sessionId}/analysis`);                const groupedFiles = analysisRes.data;
                setAnalysisData(groupedFiles);

                const fetchedTables = {};

                for (const group of Object.keys(groupedFiles)) {
                    const tsvFiles = groupedFiles[group]?.tsv || [];
                    for (const filename of tsvFiles) {

                        const tableRes = await axios.get(`/session/${sessionId}/analysis/table/${filename}`);

                        fetchedTables[filename] = tableRes.data.records || [];
                    }
                }

                setTableData(fetchedTables);
            } catch (err) {
                console.error('Failed to load analysis dashboard:', err);
                setError('Analysis data not found. Downstream analysis may have failed or is still running.');
            } finally {
                setLoading(false);
            }
        };

        fetchDashboardData();
    }, [sessionId]);

    const structured = useMemo(() => {
        if (!analysisData) return [];

        const allFiles = [];

        for (const groupKey of Object.keys(analysisData)) {
            const group = analysisData[groupKey] || {};
            const tsvFiles = group.tsv || [];
            const pngFiles = group.png || [];

            tsvFiles.forEach((filename) => {
                allFiles.push({
                    filename,
                    kind: 'tsv',
                    method: detectMethod(filename),
                    section: detectSection(filename),
                    rows: tableData[filename] || [],
                });
            });

            pngFiles.forEach((filename) => {
                allFiles.push({
                    filename,
                    kind: 'png',
                    method: detectMethod(filename),
                    section: detectSection(filename),
                });
            });
        }

        const methods = {};

        for (const file of allFiles) {
            if (!methods[file.method]) methods[file.method] = {};
            if (!methods[file.method][file.section]) methods[file.method][file.section] = [];
            methods[file.method][file.section].push(file);
        }

        return sortByOrder(Object.keys(methods), (k) => k, METHOD_ORDER).map((methodKey) => {
            const sections = methods[methodKey];
            const orderedSections = sortByOrder(Object.keys(sections), (k) => k, SECTION_ORDER).map((sectionKey) => ({
                sectionKey,
                files: sections[sectionKey].sort((a, b) => a.filename.localeCompare(b.filename)),
            }));

            return { methodKey, sections: orderedSections };
        });
    }, [analysisData, tableData]);

    if (!sessionId) return null;

    return (
        <div className="bg-white rounded-2xl shadow-md border border-gray-200 mt-6 p-6">
            <div className="flex items-start justify-between gap-4 border-b pb-4">
                <div>
                    <h3 className="text-2xl font-bold text-gray-800">Analysis Dashboard</h3>
                    <p className="text-sm text-gray-500 mt-1">
                        Structured view of summaries, metrics, tables, and figures.
                    </p>
                </div>
            </div>

            {loading && (
                <div className="mt-6 rounded-xl border border-blue-200 bg-blue-50 px-4 py-3 text-blue-700 font-medium">
                    Loading analytical data...
                </div>
            )}

            {error && (
                <div className="mt-6 rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-red-700">
                    {error}
                </div>
            )}

            {!loading && !error && structured.length > 0 && (
                <div className="mt-6 space-y-8">
                    {structured.map(({ methodKey, sections }) => (
                        <section
                            key={methodKey}
                            className="rounded-3xl border border-gray-200 bg-gray-50/70 p-5"
                        >
                            <div className="mb-5">
                                <h4 className="text-xl font-bold text-gray-800">
                                    {prettifyLabel(methodKey)}
                                </h4>
                                <p className="text-sm text-gray-500 mt-1">
                                    Organized outputs for {prettifyLabel(methodKey).toLowerCase()} analysis.
                                </p>
                            </div>

                            <div className="space-y-5">
                                {sections.map(({ sectionKey, files }) => {
                                    const pngs = files.filter((f) => f.kind === 'png');
                                    const tsvs = files.filter((f) => f.kind === 'tsv');

                                    return (
                                        <SectionCard
                                            key={`${methodKey}-${sectionKey}`}
                                            title={prettifyLabel(sectionKey)}
                                            defaultOpen={sectionKey === 'summary' || sectionKey === 'alignment_summary'}
                                            rightNote={`${tsvs.length} table${tsvs.length !== 1 ? 's' : ''}, ${pngs.length} image${pngs.length !== 1 ? 's' : ''}`}
                                        >
                                            <div className="space-y-4">
                                                {tsvs.length > 0 && (
                                                    <div className="space-y-4">
                                                        {tsvs.map((file) => (
                                                            <DataTable
                                                                key={file.filename}
                                                                title={cleanTitle(file.filename, sessionId)}
                                                                rows={file.rows}
                                                            />
                                                        ))}
                                                    </div>
                                                )}

                                                {pngs.length > 0 && (
                                                    <div className="grid grid-cols-1 xl:grid-cols-2 gap-5">
                                                        {pngs.map((file) => (
                                                            <ImageCard
                                                                key={file.filename}
                                                                sessionId={sessionId}
                                                                filename={file.filename}
                                                                sessionLabel={sessionId}
                                                            />
                                                        ))}
                                                    </div>
                                                )}
                                            </div>
                                        </SectionCard>
                                    );
                                })}
                            </div>
                        </section>
                    ))}
                </div>
            )}

            {!loading && !error && structured.length === 0 && (
                <div className="mt-6 rounded-xl border border-gray-200 bg-gray-50 px-4 py-3 text-gray-600">
                    No structured analysis outputs were found.
                </div>
            )}
        </div>
    );
}

export default AnalysisDashboard;