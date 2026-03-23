import React, { useMemo, useState } from "react";
import axios from "axios";
import { FixedSizeList as List } from "react-window";
import { useQuery } from "@tanstack/react-query";
import { motion, AnimatePresence } from "framer-motion";
import {
    Dna,
    Loader2,
    AlertCircle,
    CheckCircle2,
    XCircle,
    Minus,
    ChevronRight,
} from "lucide-react";

const CHUNK_SIZE = 80;
const ROW_HEIGHT = 118;
const METHODS = ["global", "local", "lcs"];

const api = axios.create({
    baseURL: "http://127.0.0.1:8000",
});

function parseFasta(text) {
    const lines = text
        .split("\n")
        .map((line) => line.trim())
        .filter(Boolean);

    const sequences = [];
    let current = "";

    for (const line of lines) {
        if (line.startsWith(">")) {
            if (current) sequences.push(current);
            current = "";
        } else {
            current += line;
        }
    }

    if (current) sequences.push(current);

    if (sequences.length < 2) {
        throw new Error("FASTA does not contain two valid sequences.");
    }

    return {
        querySeq: sequences[0],
        targetSeq: sequences[1],
    };
}

function buildChunks(querySeq, targetSeq, chunkSize = CHUNK_SIZE) {
    const chunks = [];

    for (let i = 0; i < querySeq.length; i += chunkSize) {
        const query = querySeq.slice(i, i + chunkSize);
        const target = targetSeq.slice(i, i + chunkSize);

        let matches = 0;
        let mismatches = 0;
        let gaps = 0;

        for (let j = 0; j < query.length; j++) {
            const q = query[j];
            const t = target[j];

            if (q === "-" || t === "-") {
                gaps += 1;
            } else if (q === t) {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }

        chunks.push({
            index: i,
            query,
            target,
            matches,
            mismatches,
            gaps,
        });
    }

    return chunks;
}

function getMatchLine(query, target) {
    return Array.from(query).map((qChar, idx) => {
        const tChar = target[idx];
        if (qChar === "-" || tChar === "-") return " ";
        if (qChar === tChar) return "|";
        return "*";
    });
}

function SequenceChar({ char, type }) {
    const baseClass =
        "inline-block min-w-[10px] text-center font-mono text-[13px] font-semibold";

    if (char === "-") {
        return <span className={`${baseClass} text-slate-300`}>-</span>;
    }

    if (type === "query") {
        return <span className={`${baseClass} text-indigo-600`}>{char}</span>;
    }

    if (type === "target") {
        return <span className={`${baseClass} text-violet-600`}>{char}</span>;
    }

    if (type === "match") {
        if (char === "|") {
            return <span className={`${baseClass} text-emerald-500`}>|</span>;
        }
        if (char === "*") {
            return <span className={`${baseClass} text-rose-400`}>*</span>;
        }
        return <span className={`${baseClass} text-transparent`}>.</span>;
    }

    return <span className={baseClass}>{char}</span>;
}

function StatPill({ icon, label, value, className = "" }) {
    return (
        <div
            className={`inline-flex items-center gap-2 rounded-full border px-3 py-1 text-xs font-medium ${className}`}
        >
            {icon}
            <span className="text-slate-600">{label}</span>
            <span className="text-slate-900">{value}</span>
        </div>
    );
}

function MethodTabs({ activeTab, onChange }) {
    return (
        <div className="inline-flex rounded-2xl border border-slate-200 bg-slate-100/80 p-1 shadow-sm">
            {METHODS.map((method) => {
                const active = activeTab === method;
                return (
                    <button
                        key={method}
                        onClick={() => onChange(method)}
                        className={`relative rounded-xl px-4 py-2 text-sm font-semibold capitalize transition-all ${
                            active
                                ? "bg-white text-slate-950 shadow-sm"
                                : "text-slate-500 hover:text-slate-800"
                        }`}
                    >
                        {active && (
                            <motion.span
                                layoutId="active-method-pill"
                                className="absolute inset-0 rounded-xl bg-white"
                                transition={{ type: "spring", stiffness: 380, damping: 30 }}
                            />
                        )}
                        <span className="relative z-10">{method}</span>
                    </button>
                );
            })}
        </div>
    );
}

function Legend() {
    return (
        <div className="flex flex-wrap items-center gap-3 text-xs">
            <StatPill
                icon={<CheckCircle2 className="h-3.5 w-3.5 text-emerald-500" />}
                label="Match"
                value="|"
                className="border-emerald-200 bg-emerald-50"
            />
            <StatPill
                icon={<XCircle className="h-3.5 w-3.5 text-rose-400" />}
                label="Mismatch"
                value="*"
                className="border-rose-200 bg-rose-50"
            />
            <StatPill
                icon={<Minus className="h-3.5 w-3.5 text-slate-400" />}
                label="Gap"
                value="-"
                className="border-slate-200 bg-slate-50"
            />
        </div>
    );
}

function AlignmentRow({ index, style, data }) {
    const chunk = data[index];
    if (!chunk) return null;

    const matchLine = getMatchLine(chunk.query, chunk.target);

    return (
        <div style={style} className="px-4 py-3">
            <div className="rounded-2xl border border-slate-200 bg-white/90 shadow-sm transition hover:shadow-md">
                <div className="flex items-center justify-between border-b border-slate-100 px-4 py-2">
                    <div className="inline-flex items-center gap-2 text-xs font-semibold tracking-wide text-slate-500 uppercase">
                        <ChevronRight className="h-3.5 w-3.5" />
                        Residues {chunk.index + 1}–{chunk.index + chunk.query.length}
                    </div>

                    <div className="flex items-center gap-2">
            <span className="rounded-full bg-emerald-50 px-2.5 py-1 text-[11px] font-semibold text-emerald-700">
              {chunk.matches} matches
            </span>
                        <span className="rounded-full bg-rose-50 px-2.5 py-1 text-[11px] font-semibold text-rose-700">
              {chunk.mismatches} mismatches
            </span>
                        <span className="rounded-full bg-slate-100 px-2.5 py-1 text-[11px] font-semibold text-slate-600">
              {chunk.gaps} gaps
            </span>
                    </div>
                </div>

                <div className="space-y-1 px-4 py-3">
                    <div className="overflow-x-auto">
                        <div className="whitespace-nowrap tracking-[0.18em]">
                            {Array.from(chunk.query).map((char, i) => (
                                <SequenceChar key={`q-${i}`} char={char} type="query" />
                            ))}
                        </div>
                    </div>

                    <div className="overflow-x-auto">
                        <div className="whitespace-nowrap tracking-[0.18em]">
                            {matchLine.map((char, i) => (
                                <SequenceChar key={`m-${i}`} char={char} type="match" />
                            ))}
                        </div>
                    </div>

                    <div className="overflow-x-auto">
                        <div className="whitespace-nowrap tracking-[0.18em]">
                            {Array.from(chunk.target).map((char, i) => (
                                <SequenceChar key={`t-${i}`} char={char} type="target" />
                            ))}
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}

export default function AlignmentViewer({ sessionId }) {
    const [activeTab, setActiveTab] = useState("global");

    const {
        data: sessionMeta,
        isLoading: metaLoading,
        isError: metaError,
    } = useQuery({
        queryKey: ["session-meta", sessionId],
        queryFn: async () => {
            const res = await api.get(`/session/${sessionId}`);
            return res.data;
        },
        enabled: Boolean(sessionId),
        staleTime: 1000 * 60 * 5,
    });

    const alignMethod = sessionMeta?.parameters?.align_method;
    const isAll = alignMethod === "all";
    const effectiveTab = isAll ? activeTab : alignMethod;

    React.useEffect(() => {
        if (!alignMethod) return;
        setActiveTab(alignMethod === "all" ? "global" : alignMethod);
    }, [alignMethod]);

    const {
        data: fastaText,
        isLoading: fileLoading,
        isError: fileError,
        error: fileErrorObj,
    } = useQuery({
        queryKey: ["alignment-file", sessionId, effectiveTab],
        queryFn: async () => {
            const filename = `${effectiveTab}_alignment.fasta`;
            const res = await api.get(`/session/${sessionId}/file/${filename}`, {
                responseType: "text",
            });
            return res.data;
        },
        enabled: Boolean(sessionId && effectiveTab),
        retry: 1,
    });

    const derived = useMemo(() => {
        if (!fastaText) {
            return {
                chunks: [],
                totalLength: 0,
                totalMatches: 0,
                totalMismatches: 0,
                totalGaps: 0,
                identity: 0,
            };
        }

        const { querySeq, targetSeq } = parseFasta(fastaText);
        const chunks = buildChunks(querySeq, targetSeq);

        const totalMatches = chunks.reduce((acc, c) => acc + c.matches, 0);
        const totalMismatches = chunks.reduce((acc, c) => acc + c.mismatches, 0);
        const totalGaps = chunks.reduce((acc, c) => acc + c.gaps, 0);
        const comparable = totalMatches + totalMismatches;
        const identity = comparable > 0 ? (totalMatches / comparable) * 100 : 0;

        return {
            chunks,
            totalLength: querySeq.length,
            totalMatches,
            totalMismatches,
            totalGaps,
            identity,
        };
    }, [fastaText]);

    if (!sessionId) return null;

    if (metaLoading) {
        return (
            <div className="mt-6 rounded-3xl border border-slate-200 bg-white p-8 shadow-sm">
                <div className="flex items-center justify-center gap-3 text-slate-500">
                    <Loader2 className="h-5 w-5 animate-spin" />
                    <span className="text-sm font-medium">Loading alignment session…</span>
                </div>
            </div>
        );
    }

    if (metaError || !sessionMeta) {
        return (
            <div className="mt-6 rounded-3xl border border-rose-200 bg-rose-50 p-6 shadow-sm">
                <div className="flex items-center gap-3 text-rose-700">
                    <AlertCircle className="h-5 w-5" />
                    <span className="font-semibold">Failed to load session metadata.</span>
                </div>
            </div>
        );
    }

    return (
        <motion.section
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="mt-6 overflow-hidden rounded-[28px] border border-slate-200 bg-white shadow-[0_10px_40px_rgba(15,23,42,0.06)]"
        >
            <div className="border-b border-slate-200 bg-gradient-to-r from-slate-50 via-white to-slate-50 px-6 py-5">
                <div className="flex flex-col gap-4 xl:flex-row xl:items-center xl:justify-between">
                    <div>
                        <div className="mb-2 inline-flex items-center gap-2 rounded-full border border-slate-200 bg-white px-3 py-1 text-xs font-semibold text-slate-600 shadow-sm">
                            <Dna className="h-3.5 w-3.5 text-indigo-500" />
                            Sequence Alignment Viewer
                        </div>

                        <h3 className="text-2xl font-bold tracking-tight text-slate-900">
                            Alignment Explorer
                        </h3>

                        <p className="mt-1 text-sm text-slate-500">
                            Inspect aligned sequence blocks with virtualization, chunk-level
                            stats, and method switching.
                        </p>
                    </div>

                    <div className="flex flex-wrap items-center gap-3">
                        {isAll && (
                            <MethodTabs activeTab={activeTab} onChange={setActiveTab} />
                        )}

                        <div className="rounded-2xl border border-slate-200 bg-white px-4 py-2 shadow-sm">
                            <div className="text-[11px] font-semibold uppercase tracking-wide text-slate-400">
                                Method
                            </div>
                            <div className="text-sm font-bold capitalize text-slate-900">
                                {effectiveTab}
                            </div>
                        </div>
                    </div>
                </div>

                <div className="mt-4 flex flex-wrap items-center gap-3">
                    <StatPill
                        icon={<Dna className="h-3.5 w-3.5 text-slate-500" />}
                        label="Length"
                        value={`${derived.totalLength} bp`}
                        className="bg-slate-50"
                    />
                    <StatPill
                        icon={<CheckCircle2 className="h-3.5 w-3.5 text-emerald-500" />}
                        label="Matches"
                        value={derived.totalMatches}
                        className="bg-emerald-50"
                    />
                    <StatPill
                        icon={<XCircle className="h-3.5 w-3.5 text-rose-400" />}
                        label="Mismatches"
                        value={derived.totalMismatches}
                        className="bg-rose-50"
                    />
                    <StatPill
                        icon={<Minus className="h-3.5 w-3.5 text-slate-400" />}
                        label="Gaps"
                        value={derived.totalGaps}
                        className="bg-slate-50"
                    />
                    <StatPill
                        icon={<CheckCircle2 className="h-3.5 w-3.5 text-indigo-500" />}
                        label="Identity"
                        value={`${derived.identity.toFixed(2)}%`}
                        className="bg-indigo-50"
                    />
                </div>

                <div className="mt-4">
                    <Legend />
                </div>
            </div>

            <div className="bg-slate-50/60 p-4">
                <AnimatePresence mode="wait">
                    {fileLoading ? (
                        <motion.div
                            key="loading"
                            initial={{ opacity: 0 }}
                            animate={{ opacity: 1 }}
                            exit={{ opacity: 0 }}
                            className="flex h-[420px] items-center justify-center rounded-3xl border border-dashed border-slate-300 bg-white"
                        >
                            <div className="flex items-center gap-3 text-slate-500">
                                <Loader2 className="h-5 w-5 animate-spin" />
                                <span className="text-sm font-medium">
                  Parsing FASTA alignment…
                </span>
                            </div>
                        </motion.div>
                    ) : fileError ? (
                        <motion.div
                            key="error"
                            initial={{ opacity: 0 }}
                            animate={{ opacity: 1 }}
                            exit={{ opacity: 0 }}
                            className="flex h-[420px] items-center justify-center rounded-3xl border border-rose-200 bg-rose-50 p-6"
                        >
                            <div className="max-w-md text-center">
                                <AlertCircle className="mx-auto mb-3 h-8 w-8 text-rose-500" />
                                <h4 className="text-base font-bold text-rose-700">
                                    Could not load alignment data
                                </h4>
                                <p className="mt-2 text-sm text-rose-600">
                                    {fileErrorObj?.message ||
                                        `The ${effectiveTab} alignment file may still be generating.`}
                                </p>
                            </div>
                        </motion.div>
                    ) : derived.chunks.length > 0 ? (
                        <motion.div
                            key={effectiveTab}
                            initial={{ opacity: 0, y: 8 }}
                            animate={{ opacity: 1, y: 0 }}
                            exit={{ opacity: 0, y: -8 }}
                            className="overflow-hidden rounded-3xl border border-slate-200 bg-white shadow-sm"
                        >
                            <List
                                height={420}
                                itemCount={derived.chunks.length}
                                itemSize={ROW_HEIGHT}
                                width="100%"
                                itemData={derived.chunks}
                            >
                                {AlignmentRow}
                            </List>
                        </motion.div>
                    ) : (
                        <motion.div
                            key="empty"
                            initial={{ opacity: 0 }}
                            animate={{ opacity: 1 }}
                            exit={{ opacity: 0 }}
                            className="flex h-[420px] items-center justify-center rounded-3xl border border-slate-200 bg-white"
                        >
                            <div className="text-center text-slate-500">
                                <Dna className="mx-auto mb-3 h-8 w-8 text-slate-400" />
                                <p className="font-medium">No alignment chunks available.</p>
                            </div>
                        </motion.div>
                    )}
                </AnimatePresence>
            </div>
        </motion.section>
    );
}