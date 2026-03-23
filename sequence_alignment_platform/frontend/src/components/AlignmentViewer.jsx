import React, { useEffect, useMemo, useState } from "react";
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
} from "lucide-react";

const CHUNK_SIZE = 60;
const ROW_HEIGHT = 108;
const METHODS = ["global", "local", "lcs"];

const api = axios.create({
    baseURL: "http://127.0.0.1:8000",
});

function parseFasta(text) {
    const lines = text
        .split("\n")
        .map((line) => line.trim())
        .filter(Boolean);

    const entries = [];
    let currentHeader = "";
    let currentSequence = "";

    for (const line of lines) {
        if (line.startsWith(">")) {
            if (currentSequence) {
                entries.push({
                    header: currentHeader,
                    sequence: currentSequence,
                });
            }
            currentHeader = line.slice(1).trim();
            currentSequence = "";
        } else {
            currentSequence += line;
        }
    }

    if (currentSequence) {
        entries.push({
            header: currentHeader,
            sequence: currentSequence,
        });
    }

    if (entries.length < 2) {
        throw new Error("FASTA does not contain two valid sequences.");
    }

    return {
        queryHeader: entries[0].header || "Query",
        targetHeader: entries[1].header || "Target",
        querySeq: entries[0].sequence,
        targetSeq: entries[1].sequence,
    };
}

function getMatchSymbol(queryChar, targetChar) {
    if (queryChar === "-" || targetChar === "-") return " ";
    if (queryChar === targetChar) return "|";
    return ".";
}

function nonGapCount(sequence) {
    return (sequence.match(/[^-]/g) || []).length;
}

function padLeft(value, width) {
    return String(value).padStart(width, " ");
}

function padRight(value, width) {
    return String(value).padEnd(width, " ");
}

function buildChunks(querySeq, targetSeq, chunkSize = CHUNK_SIZE) {
    const chunks = [];
    let queryPosition = 1;
    let targetPosition = 1;

    for (let i = 0; i < querySeq.length; i += chunkSize) {
        const query = querySeq.slice(i, i + chunkSize);
        const target = targetSeq.slice(i, i + chunkSize);

        let matchLine = "";
        let matches = 0;
        let mismatches = 0;
        let gaps = 0;

        for (let j = 0; j < query.length; j++) {
            const symbol = getMatchSymbol(query[j], target[j]);
            matchLine += symbol;

            if (symbol === "|") matches += 1;
            else if (symbol === ".") mismatches += 1;
            else gaps += 1;
        }

        const queryStart = queryPosition;
        const targetStart = targetPosition;

        queryPosition += nonGapCount(query);
        targetPosition += nonGapCount(target);

        chunks.push({
            query,
            target,
            matchLine,
            queryStart,
            queryEnd: queryPosition - 1,
            targetStart,
            targetEnd: targetPosition - 1,
            matches,
            mismatches,
            gaps,
        });
    }

    return chunks;
}

function StatPill({ label, value, className = "" }) {
    return (
        <div className={`inline-flex items-center rounded-full border px-3 py-1 text-xs font-medium ${className}`}>
            <span className="text-slate-600">{label}</span>
            <span className="ml-2 text-slate-900">{value}</span>
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
                        type="button"
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

function AlignmentRow({ index, style, data }) {
    const {
        chunks,
        queryLabel,
        targetLabel,
    } = data;

    const chunk = chunks[index];
    if (!chunk) return null;

    const queryLine = `${padRight(queryLabel, 12)} ${padLeft(chunk.queryStart, 6)}  ${chunk.query}  ${padLeft(chunk.queryEnd, 6)}`;
    const matchLine = `${padRight("", 12)} ${padLeft("", 6)}  ${chunk.matchLine}`;
    const targetLine = `${padRight(targetLabel, 12)} ${padLeft(chunk.targetStart, 6)}  ${chunk.target}  ${padLeft(chunk.targetEnd, 6)}`;

    return (
        <div style={style} className="px-4 py-2">
            <div className="rounded-2xl border border-slate-200 bg-white px-4 py-3 shadow-sm">
        <pre
            className="m-0 overflow-x-auto text-[13px] leading-6 text-slate-800"
            style={{
                fontFamily:
                    'ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace',
                fontVariantLigatures: "none",
                whiteSpace: "pre",
            }}
        >
          <span className="text-indigo-700">{queryLine}</span>
            {"\n"}
            <span className="text-emerald-600">{matchLine}</span>
            {"\n"}
            <span className="text-violet-700">{targetLine}</span>
        </pre>
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
            const response = await api.get(`/session/${sessionId}`);
            return response.data;
        },
        enabled: Boolean(sessionId),
        staleTime: 5 * 60 * 1000,
    });

    const alignMethod = sessionMeta?.parameters?.align_method;
    const isAllMethods = alignMethod === "all";
    const effectiveTab = isAllMethods ? activeTab : alignMethod;

    useEffect(() => {
        if (!alignMethod) return;
        setActiveTab(alignMethod === "all" ? "global" : alignMethod);
    }, [alignMethod]);

    const {
        data: fastaText,
        isLoading: fileLoading,
        isError: fileError,
        error: fileErrorObject,
    } = useQuery({
        queryKey: ["alignment-file", sessionId, effectiveTab],
        queryFn: async () => {
            const filename = `${effectiveTab}_alignment.fasta`;
            const response = await api.get(`/session/${sessionId}/file/${filename}`, {
                responseType: "text",
            });
            return response.data;
        },
        enabled: Boolean(sessionId && effectiveTab),
        retry: 1,
    });

    const derived = useMemo(() => {
        if (!fastaText) {
            return {
                queryHeader: "Query",
                targetHeader: "Target",
                chunks: [],
                totalLength: 0,
                totalMatches: 0,
                totalMismatches: 0,
                totalGaps: 0,
                identity: 0,
            };
        }

        const { queryHeader, targetHeader, querySeq, targetSeq } = parseFasta(fastaText);
        const chunks = buildChunks(querySeq, targetSeq, CHUNK_SIZE);

        const totalMatches = chunks.reduce((sum, chunk) => sum + chunk.matches, 0);
        const totalMismatches = chunks.reduce((sum, chunk) => sum + chunk.mismatches, 0);
        const totalGaps = chunks.reduce((sum, chunk) => sum + chunk.gaps, 0);
        const comparableSites = totalMatches + totalMismatches;
        const identity = comparableSites > 0 ? (totalMatches / comparableSites) * 100 : 0;

        return {
            queryHeader,
            targetHeader,
            chunks,
            totalLength: Math.max(querySeq.length, targetSeq.length),
            totalMatches,
            totalMismatches,
            totalGaps,
            identity,
        };
    }, [fastaText]);

    const listData = useMemo(
        () => ({
            chunks: derived.chunks,
            queryLabel: "Query",
            targetLabel: "Target",
        }),
        [derived.chunks]
    );

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
                        <h3 className="text-2xl font-bold tracking-tight text-slate-900">
                            Alignment Explorer
                        </h3>
                    </div>

                    <div className="flex flex-wrap items-center gap-3">
                        {isAllMethods && (
                            <MethodTabs activeTab={activeTab} onChange={setActiveTab} />
                        )}

                        <div className="rounded-2xl border border-slate-200 bg-white px-4 py-2 shadow-sm">
                        </div>
                    </div>
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
                                <span className="text-sm font-medium">Parsing FASTA alignment…</span>
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
                                    {fileErrorObject?.message ||
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
                                itemData={listData}
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