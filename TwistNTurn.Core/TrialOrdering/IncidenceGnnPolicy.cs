using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Threading;

namespace TwistNTurn.TrialOrdering
{
    public sealed class LinearWeights
    {
        public readonly int OutDim;
        public readonly int InDim;
        public readonly float[] W;  // [OutDim * InDim] row-major
        public readonly float[] WT; // [InDim * OutDim] transposed for SIMD register-accumulated GEMV
        public readonly float[] B;  // [OutDim]

        public LinearWeights(BinaryReader br)
        {
            OutDim = br.ReadInt32();
            InDim = br.ReadInt32();
            W = new float[OutDim * InDim];
            WT = new float[InDim * OutDim];
            B = new float[OutDim];
            byte[] wBytes = br.ReadBytes(W.Length * sizeof(float));
            if (wBytes.Length != W.Length * sizeof(float))
                throw new EndOfStreamException("Unexpected end of LDLM weight tensor.");
            Buffer.BlockCopy(wBytes, 0, W, 0, wBytes.Length);
            byte[] bBytes = br.ReadBytes(B.Length * sizeof(float));
            if (bBytes.Length != B.Length * sizeof(float))
                throw new EndOfStreamException("Unexpected end of LDLM bias tensor.");
            Buffer.BlockCopy(bBytes, 0, B, 0, bBytes.Length);

            for (int o = 0; o < OutDim; o++)
                for (int i = 0; i < InDim; i++)
                    WT[i * OutDim + o] = W[o * InDim + i];
        }
    }

    /// <summary>
    /// Immutable, thread-safe weights for a TwistNTurn IncidenceGNN (LDLM v1/v2).
    /// Load once via <see cref="Default"/> (embedded resource) and share across policies.
    /// </summary>
    public sealed class IncidenceGnnModel
    {
        public const uint Magic = 0x4C444C4D; // "LDLM"
        public const string DefaultResourceName = "TwistNTurn.TrialOrdering.Weights.tnt_side2_h16_l1.ldlm";

        private static readonly Lazy<IncidenceGnnModel> s_default =
            new Lazy<IncidenceGnnModel>(LoadDefaultEmbedded, LazyThreadSafetyMode.ExecutionAndPublication);

        /// <summary>
        /// Thread-safe singleton for the embedded 21.3 KB unified side-predicting model (tnt_side2_h16_l1).
        /// </summary>
        public static IncidenceGnnModel Default => s_default.Value;

        public readonly int Version;
        public readonly int Fe, Fc, Fv, Fg;
        public readonly int Hidden;
        public readonly int NumLayers;
        public readonly bool UseColorPool;
        public readonly bool UseMaxPool;

        public readonly LinearWeights EIn, CIn, VIn;
        public readonly LinearWeights[] VUpd, CUpd, EUpd;
        public readonly float[]?[] VNorm, CNorm, ENorm;
        public readonly LinearWeights EHead0, EHead1;
        public readonly LinearWeights CHead0, CHead1;
        public bool HasSideHead => EHead1.OutDim >= 2 && CHead1.OutDim >= 2;

        public IncidenceGnnModel(string path)
            : this(OpenRead(path), leaveOpen: false)
        {
        }

        public IncidenceGnnModel(Stream stream, bool leaveOpen = true)
        {
            ArgumentNullException.ThrowIfNull(stream);
            using BinaryReader br = new BinaryReader(stream, System.Text.Encoding.UTF8, leaveOpen);
            uint magic = br.ReadUInt32();
            Version = br.ReadInt32();
            if (magic != Magic || (Version != 1 && Version != 2))
                throw new InvalidDataException($"Invalid LDLM model stream: magic=0x{magic:X8}, version={Version}");

            Fe = br.ReadInt32();
            Fc = br.ReadInt32();
            Fv = br.ReadInt32();
            Fg = br.ReadInt32();
            if (Fe != GraphFeatureExtractor.EdgeFeatDim || Fc != GraphFeatureExtractor.CellFeatDim ||
                Fv != GraphFeatureExtractor.VertFeatDim || Fg != GraphFeatureExtractor.GlobalFeatDim)
            {
                throw new InvalidDataException($"LDLM feature dimensions ({Fe},{Fc},{Fv},{Fg}) do not match extractor ({GraphFeatureExtractor.EdgeFeatDim},{GraphFeatureExtractor.CellFeatDim},{GraphFeatureExtractor.VertFeatDim},{GraphFeatureExtractor.GlobalFeatDim}).");
            }

            Hidden = br.ReadInt32();
            NumLayers = br.ReadInt32();
            UseColorPool = br.ReadInt32() != 0;
            UseMaxPool = Version >= 2 && br.ReadInt32() != 0;

            EIn = new LinearWeights(br);
            CIn = new LinearWeights(br);
            VIn = new LinearWeights(br);

            VUpd = new LinearWeights[NumLayers];
            CUpd = new LinearWeights[NumLayers];
            EUpd = new LinearWeights[NumLayers];
            VNorm = new float[NumLayers][];
            CNorm = new float[NumLayers][];
            ENorm = new float[NumLayers][];
            for (int l = 0; l < NumLayers; l++)
            {
                VUpd[l] = new LinearWeights(br);
                if (Version >= 2) VNorm[l] = ReadFloatArray(br, Hidden);
                CUpd[l] = new LinearWeights(br);
                if (Version >= 2) CNorm[l] = ReadFloatArray(br, Hidden);
                EUpd[l] = new LinearWeights(br);
                if (Version >= 2) ENorm[l] = ReadFloatArray(br, Hidden);
            }

            EHead0 = new LinearWeights(br);
            EHead1 = new LinearWeights(br);
            CHead0 = new LinearWeights(br);
            CHead1 = new LinearWeights(br);
        }

        private static FileStream OpenRead(string path)
        {
            return new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read);
        }

        private static IncidenceGnnModel LoadDefaultEmbedded()
        {
            using Stream? stream = typeof(IncidenceGnnModel).Assembly.GetManifestResourceStream(DefaultResourceName);
            if (stream == null)
                throw new InvalidOperationException($"Embedded LDLM model resource '{DefaultResourceName}' was not found.");
            return new IncidenceGnnModel(stream, leaveOpen: true);
        }

        private static float[] ReadFloatArray(BinaryReader br, int n)
        {
            float[] arr = new float[n];
            byte[] bytes = br.ReadBytes(n * sizeof(float));
            if (bytes.Length != n * sizeof(float))
                throw new EndOfStreamException("Unexpected end of LDLM norm weights.");
            Buffer.BlockCopy(bytes, 0, arr, 0, bytes.Length);
            return arr;
        }
    }

    /// <summary>
    /// Zero-allocation pure-C# SIMD inference engine and <see cref="ITrialOrderingPolicy"/>
    /// for TwistNTurn. Fully AOT-compatible for .NET desktop and Blazor WASM AOT.
    /// Instances hold mutable per-solve scratch buffers and must not be shared concurrently across threads.
    /// </summary>
    public sealed class IncidenceGnnPolicy : ITrialOrderingPolicy
    {
        private readonly IncidenceGnnModel _m;
        private readonly PostProgressAction _postAction;
        private readonly int _reorderThreshold;
        private readonly int _minActiveTrialsToScore;
        private readonly bool _useSmartCooldown;
        private readonly int _rescoreEdgeDelta;
        private readonly bool _autoMode;
        private readonly bool _useSidePrediction;

        private Mesh? _cachedMesh;
        private MeshTopologyCache? _topo;
        private int _lastFullNnEmptyEdges;
        private int _lastRecDepth = -1;
        private int _lastConfigHash;
        private bool _hasCachedScores;

        private byte[] _lastEdgeState = Array.Empty<byte>();
        private bool[] _edgeFailed = Array.Empty<bool>();
        private bool[] _cellFailed = Array.Empty<bool>();
        private bool[] _edgeDirty = Array.Empty<bool>();
        private bool[] _edgeDirty2 = Array.Empty<bool>();
        private bool _hasSnapshotState;

        private readonly float[] _gFeats = new float[GraphFeatureExtractor.GlobalFeatDim];
        private float[] _eFeats = Array.Empty<float>();
        private float[] _cFeats = Array.Empty<float>();
        private float[] _vFeats = Array.Empty<float>();
        private int[] _eColors = Array.Empty<int>();
        private int[] _cColors = Array.Empty<int>();
        private int[] _vEsScratch = Array.Empty<int>();

        private float[] _hE = Array.Empty<float>();
        private float[] _hC = Array.Empty<float>();
        private float[] _hV = Array.Empty<float>();
        private float[] _hENext = Array.Empty<float>();
        private float[] _hCNext = Array.Empty<float>();
        private float[] _hVNext = Array.Empty<float>();

        private float[] _vMsg = Array.Empty<float>();
        private float[] _vMax = Array.Empty<float>();
        private float[] _cMsg = Array.Empty<float>();
        private float[] _cMax = Array.Empty<float>();
        private float[] _vCnt = Array.Empty<float>();
        private float[] _cCnt = Array.Empty<float>();

        private float[] _posSum = Array.Empty<float>();
        private float[] _negSum = Array.Empty<float>();
        private float[] _posCnt = Array.Empty<float>();
        private float[] _negCnt = Array.Empty<float>();

        private float[] _eScores = Array.Empty<float>();
        private float[] _cScores = Array.Empty<float>();
        private float[] _eSides = Array.Empty<float>();
        private float[] _cSides = Array.Empty<float>();
        private float[] _concatTmp = Array.Empty<float>();

        private float[] _sortScores = Array.Empty<float>();
        private int[] _sortIndices = Array.Empty<int>();
        private IAction[] _sortScratch = Array.Empty<IAction>();

        public IncidenceGnnPolicy(
            IncidenceGnnModel? model = null,
            PostProgressAction postAction = PostProgressAction.Continue,
            int reorderThreshold = 1,
            int minActiveTrialsToScore = 8,
            bool useSmartCooldown = false,
            int rescoreEdgeDelta = 0,
            bool autoMode = true,
            bool useSidePrediction = true)
        {
            _m = model ?? IncidenceGnnModel.Default;
            _postAction = postAction;
            _reorderThreshold = reorderThreshold;
            _minActiveTrialsToScore = minActiveTrialsToScore;
            _useSmartCooldown = useSmartCooldown || autoMode;
            _rescoreEdgeDelta = rescoreEdgeDelta;
            _autoMode = autoMode;
            _useSidePrediction = useSidePrediction;
        }

        /// <summary>
        /// Creates a new per-solve policy backed by the embedded default model in automatic mode.
        /// </summary>
        public static IncidenceGnnPolicy CreateDefault() => new IncidenceGnnPolicy(IncidenceGnnModel.Default, autoMode: true);

        public float[] EdgeScores => _eScores;
        public float[] CellScores => _cScores;
        public float[] EdgeSides => _eSides;
        public float[] CellSides => _cSides;

        public void OnSolveStart(Mesh mesh)
        {
            _hasCachedScores = false;
            _hasSnapshotState = false;
            _lastRecDepth = -1;
        }

        private static int ComputeConfigHash(Mesh mesh)
        {
            int h = mesh.IterativeSolverDepth;
            h = (h * 397) ^ mesh.IterativeRecMaxDepth;
            int flags = (mesh.ConsiderMultipleLoops ? 1 : 0)
                      | (mesh.UseColoring ? 2 : 0)
                      | (mesh.UseCellColoring ? 4 : 0)
                      | (mesh.UseEdgeRestricts ? 8 : 0)
                      | (mesh.UseMerging ? 16 : 0)
                      | ((mesh.UseCellPairs || mesh.UseCellPairsTopLevel) ? 32 : 0);
            return (h * 397) ^ flags;
        }

        private void EnsureBuffers(Mesh mesh)
        {
            if (!ReferenceEquals(_cachedMesh, mesh) || _topo == null
                || _topo.NumEdges != mesh.Edges.Count || _topo.NumCells != mesh.Cells.Count || _topo.NumVerts != mesh.Intersections.Count)
            {
                _cachedMesh = mesh;
                _topo = new MeshTopologyCache(mesh);
                _hasSnapshotState = false;
                _hasCachedScores = false;
                _lastRecDepth = -1;
            }

            int V = _topo.NumVerts, E = _topo.NumEdges, C = _topo.NumCells, H = _m.Hidden;
            if (_lastEdgeState.Length < E)
            {
                _lastEdgeState = new byte[E];
                _edgeFailed = new bool[E];
                _edgeDirty = new bool[E];
                _edgeDirty2 = new bool[E];
            }
            if (_cellFailed.Length < C) _cellFailed = new bool[C];
            if (_eFeats.Length < E * GraphFeatureExtractor.EdgeFeatDim) _eFeats = new float[E * GraphFeatureExtractor.EdgeFeatDim];
            if (_cFeats.Length < C * GraphFeatureExtractor.CellFeatDim) _cFeats = new float[C * GraphFeatureExtractor.CellFeatDim];
            if (_vFeats.Length < V * GraphFeatureExtractor.VertFeatDim) _vFeats = new float[V * GraphFeatureExtractor.VertFeatDim];
            if (_eColors.Length < E) _eColors = new int[E];
            if (_cColors.Length < C) _cColors = new int[C];
            if (_vEsScratch.Length < V) _vEsScratch = new int[V];

            if (_hE.Length < E * H) { _hE = new float[E * H]; _hENext = new float[E * H]; }
            if (_hC.Length < C * H) { _hC = new float[C * H]; _hCNext = new float[C * H]; }
            if (_hV.Length < V * H) { _hV = new float[V * H]; _hVNext = new float[V * H]; }
            if (_vMsg.Length < V * H) { _vMsg = new float[V * H]; _vMax = new float[V * H]; _vCnt = new float[V]; }
            if (_cMsg.Length < C * H) { _cMsg = new float[C * H]; _cMax = new float[C * H]; _cCnt = new float[C]; }

            if (_eScores.Length < E) { _eScores = new float[E]; _eSides = new float[E]; }
            if (_cScores.Length < C) { _cScores = new float[C]; _cSides = new float[C]; }
            if (_concatTmp.Length < 5 * H + 32) _concatTmp = new float[5 * H + 32];
        }

        private void EnsureColorPoolBuffers(int numGroups, int H)
        {
            if (_posSum.Length < numGroups * H)
            {
                int cap = Math.Max(numGroups * 2, 32);
                _posSum = new float[cap * H];
                _negSum = new float[cap * H];
                _posCnt = new float[cap];
                _negCnt = new float[cap];
            }
        }

        public void ScoreAllNodes(Mesh mesh, int iterativeRecDepth)
        {
            EnsureBuffers(mesh);
            MeshTopologyCache topo = _topo!;
            int V = topo.NumVerts, E = topo.NumEdges, C = topo.NumCells, H = _m.Hidden;

            GraphFeatureExtractor.Extract(mesh, topo, _gFeats, _eFeats, _cFeats, _vFeats, _eColors, _cColors, _vEsScratch, iterativeRecDepth);

            ProjectInputs(_m.EIn, _eFeats, _m.Fe, _gFeats, _m.Fg, E, _hE);
            ProjectInputs(_m.CIn, _cFeats, _m.Fc, _gFeats, _m.Fg, C, _hC);
            ProjectInputs(_m.VIn, _vFeats, _m.Fv, _gFeats, _m.Fg, V, _hV);

            RunMessagePassing(topo, V, E, C, H);

            EvalEdgeHead(_m.EHead0, _m.EHead1, _hE, mesh.Edges, E, H, _eScores, _eSides);
            EvalCellHead(_m.CHead0, _m.CHead1, _hC, mesh.Cells, C, H, _cScores, _cSides);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void AddVec(float[] dst, int dOff, float[] src, int sOff, int H)
        {
            ref float dRef = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(dst), dOff);
            ref float sRef = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(src), sOff);
            int h = 0;
            if (Vector128.IsHardwareAccelerated)
            {
                for (; h <= H - 4; h += 4)
                    (Vector128.LoadUnsafe(ref dRef, (nuint)h) + Vector128.LoadUnsafe(ref sRef, (nuint)h)).StoreUnsafe(ref dRef, (nuint)h);
            }
            for (; h < H; h++)
                Unsafe.Add(ref dRef, h) += Unsafe.Add(ref sRef, h);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void AddAndMaxVec(float[] sumDst, float[] maxDst, int dOff, float[] src, int sOff, int H)
        {
            ref float sumRef = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(sumDst), dOff);
            ref float maxRef = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(maxDst), dOff);
            ref float sRef = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(src), sOff);
            int h = 0;
            if (Vector128.IsHardwareAccelerated)
            {
                for (; h <= H - 4; h += 4)
                {
                    Vector128<float> sv = Vector128.LoadUnsafe(ref sRef, (nuint)h);
                    (Vector128.LoadUnsafe(ref sumRef, (nuint)h) + sv).StoreUnsafe(ref sumRef, (nuint)h);
                    Vector128.Max(Vector128.LoadUnsafe(ref maxRef, (nuint)h), sv).StoreUnsafe(ref maxRef, (nuint)h);
                }
            }
            for (; h < H; h++)
            {
                float v = Unsafe.Add(ref sRef, h);
                Unsafe.Add(ref sumRef, h) += v;
                if (v > Unsafe.Add(ref maxRef, h)) Unsafe.Add(ref maxRef, h) = v;
            }
        }

        private void RunMessagePassing(MeshTopologyCache topo, int V, int E, int C, int H)
        {
            bool useMax = _m.UseMaxPool;
            for (int l = 0; l < _m.NumLayers; l++)
            {
                Array.Clear(_vMsg, 0, V * H);
                Array.Clear(_vCnt, 0, V);
                if (useMax) Array.Fill(_vMax, -1e9f, 0, V * H);
                for (int e = 0; e < E; e++)
                {
                    int v0 = topo.EdgeV0[e], v1 = topo.EdgeV1[e];
                    int eo = e * H, vo0 = v0 * H, vo1 = v1 * H;
                    if (useMax)
                    {
                        AddAndMaxVec(_vMsg, _vMax, vo0, _hE, eo, H);
                        AddAndMaxVec(_vMsg, _vMax, vo1, _hE, eo, H);
                    }
                    else
                    {
                        AddVec(_vMsg, vo0, _hE, eo, H);
                        AddVec(_vMsg, vo1, _hE, eo, H);
                    }
                    _vCnt[v0] += 1f;
                    _vCnt[v1] += 1f;
                }
                for (int v = 0; v < V; v++)
                {
                    float inv = _vCnt[v] > 0f ? 1f / _vCnt[v] : 0f;
                    int vo = v * H;
                    for (int h = 0; h < H; h++)
                    {
                        _concatTmp[h] = _hV[vo + h];
                        _concatTmp[H + h] = _vMsg[vo + h] * inv;
                        if (useMax) _concatTmp[2 * H + h] = _vMax[vo + h] > -1e8f ? _vMax[vo + h] : 0f;
                    }
                    LinearAddResidualRelu(_m.VUpd[l], _concatTmp, _hV, vo, _hVNext, vo);
                    if (_m.VNorm[l] != null) ApplyRmsNorm(_hVNext, vo, H, _m.VNorm[l]!);
                }
                (_hV, _hVNext) = (_hVNext, _hV);

                Array.Clear(_cMsg, 0, C * H);
                Array.Clear(_cCnt, 0, C);
                if (useMax) Array.Fill(_cMax, -1e9f, 0, C * H);
                for (int e = 0; e < E; e++)
                {
                    int c0 = topo.EdgeC0[e], c1 = topo.EdgeC1[e];
                    int eo = e * H;
                    if (c0 >= 0)
                    {
                        if (useMax) AddAndMaxVec(_cMsg, _cMax, c0 * H, _hE, eo, H);
                        else AddVec(_cMsg, c0 * H, _hE, eo, H);
                        _cCnt[c0] += 1f;
                    }
                    if (c1 >= 0)
                    {
                        if (useMax) AddAndMaxVec(_cMsg, _cMax, c1 * H, _hE, eo, H);
                        else AddVec(_cMsg, c1 * H, _hE, eo, H);
                        _cCnt[c1] += 1f;
                    }
                }

                int maxCCol = _m.UseColorPool ? BuildColorPool(_hC, _cColors, C, H) : 0;
                int cColOff = (useMax ? 3 : 2) * H;
                for (int c = 0; c < C; c++)
                {
                    float inv = _cCnt[c] > 0f ? 1f / _cCnt[c] : 0f;
                    int co = c * H;
                    for (int h = 0; h < H; h++)
                    {
                        _concatTmp[h] = _hC[co + h];
                        _concatTmp[H + h] = _cMsg[co + h] * inv;
                        if (useMax) _concatTmp[2 * H + h] = _cMax[co + h] > -1e8f ? _cMax[co + h] : 0f;
                    }
                    if (_m.UseColorPool) FillColorPoolMsg(_cColors[c], maxCCol, H, _concatTmp, cColOff);
                    LinearAddResidualRelu(_m.CUpd[l], _concatTmp, _hC, co, _hCNext, co);
                    if (_m.CNorm[l] != null) ApplyRmsNorm(_hCNext, co, H, _m.CNorm[l]!);
                }
                (_hC, _hCNext) = (_hCNext, _hC);

                int maxECol = _m.UseColorPool ? BuildColorPool(_hE, _eColors, E, H) : 0;
                for (int e = 0; e < E; e++)
                {
                    int eo = e * H;
                    int vo0 = topo.EdgeV0[e] * H;
                    int vo1 = topo.EdgeV1[e] * H;
                    int c0 = topo.EdgeC0[e];
                    int c1 = topo.EdgeC1[e];
                    for (int h = 0; h < H; h++)
                    {
                        _concatTmp[h] = _hE[eo + h];
                        _concatTmp[H + h] = 0.5f * (_hV[vo0 + h] + _hV[vo1 + h]);
                        float cVal = 0f;
                        if (c0 >= 0 && c1 >= 0) cVal = 0.5f * (_hC[c0 * H + h] + _hC[c1 * H + h]);
                        else if (c0 >= 0) cVal = _hC[c0 * H + h];
                        else if (c1 >= 0) cVal = _hC[c1 * H + h];
                        _concatTmp[2 * H + h] = cVal;
                    }
                    if (_m.UseColorPool) FillColorPoolMsg(_eColors[e], maxECol, H, _concatTmp, 3 * H);
                    LinearAddResidualRelu(_m.EUpd[l], _concatTmp, _hE, eo, _hENext, eo);
                    if (_m.ENorm[l] != null) ApplyRmsNorm(_hENext, eo, H, _m.ENorm[l]!);
                }
                (_hE, _hENext) = (_hENext, _hE);
            }
        }

        private int BuildColorPool(float[] hArr, int[] colors, int N, int H)
        {
            int maxCol = 0;
            for (int i = 0; i < N; i++) { int a = Math.Abs(colors[i]); if (a > maxCol) maxCol = a; }
            if (maxCol == 0) return 0;

            int numGroups = maxCol + 1;
            EnsureColorPoolBuffers(numGroups, H);
            Array.Clear(_posSum, 0, numGroups * H);
            Array.Clear(_negSum, 0, numGroups * H);
            Array.Clear(_posCnt, 0, numGroups);
            Array.Clear(_negCnt, 0, numGroups);

            for (int i = 0; i < N; i++)
            {
                int col = colors[i];
                if (col > 0)
                {
                    int go = col * H, io = i * H;
                    for (int h = 0; h < H; h++) _posSum[go + h] += hArr[io + h];
                    _posCnt[col] += 1f;
                }
                else if (col < 0)
                {
                    int cid = -col;
                    int go = cid * H, io = i * H;
                    for (int h = 0; h < H; h++) _negSum[go + h] += hArr[io + h];
                    _negCnt[cid] += 1f;
                }
            }
            for (int g = 1; g < numGroups; g++)
            {
                int go = g * H;
                if (_posCnt[g] > 1f) { float inv = 1f / _posCnt[g]; for (int h = 0; h < H; h++) _posSum[go + h] *= inv; }
                if (_negCnt[g] > 1f) { float inv = 1f / _negCnt[g]; for (int h = 0; h < H; h++) _negSum[go + h] *= inv; }
            }
            return maxCol;
        }

        private void FillColorPoolMsg(int col, int maxCol, int H, float[] dst, int dstOff)
        {
            if (col == 0 || maxCol == 0) Array.Clear(dst, dstOff, 2 * H);
            else if (col > 0) { int go = col * H; Array.Copy(_posSum, go, dst, dstOff, H); Array.Copy(_negSum, go, dst, dstOff + H, H); }
            else { int go = (-col) * H; Array.Copy(_negSum, go, dst, dstOff, H); Array.Copy(_posSum, go, dst, dstOff + H, H); }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ApplyRmsNorm(float[] arr, int off, int H, float[] gamma)
        {
            float sumSq = 0f;
            for (int h = 0; h < H; h++) { float v = arr[off + h]; sumSq += v * v; }
            float invRms = 1f / MathF.Sqrt(sumSq / H + 1e-5f);
            for (int h = 0; h < H; h++) arr[off + h] = arr[off + h] * invRms * gamma[h];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void AccumulateMatVec(float[] wt, int wtRowStart, float[] x, int xOff, int len, Span<float> acc, int outDim)
        {
            ref float wBase = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(wt), wtRowStart * outDim);
            ref float xBase = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(x), xOff);
            ref float aRef = ref MemoryMarshal.GetReference(acc);

            if (Vector128.IsHardwareAccelerated && outDim == 16)
            {
                Vector128<float> a0 = Vector128.LoadUnsafe(ref aRef, 0); Vector128<float> a1 = Vector128.LoadUnsafe(ref aRef, 4);
                Vector128<float> a2 = Vector128.LoadUnsafe(ref aRef, 8); Vector128<float> a3 = Vector128.LoadUnsafe(ref aRef, 12);
                for (int j = 0; j < len; j++)
                {
                    Vector128<float> vx = Vector128.Create(Unsafe.Add(ref xBase, j));
                    ref float wRow = ref Unsafe.Add(ref wBase, j * 16);
                    a0 += vx * Vector128.LoadUnsafe(ref wRow, 0); a1 += vx * Vector128.LoadUnsafe(ref wRow, 4);
                    a2 += vx * Vector128.LoadUnsafe(ref wRow, 8); a3 += vx * Vector128.LoadUnsafe(ref wRow, 12);
                }
                a0.StoreUnsafe(ref aRef, 0); a1.StoreUnsafe(ref aRef, 4); a2.StoreUnsafe(ref aRef, 8); a3.StoreUnsafe(ref aRef, 12);
                return;
            }
            if (Vector128.IsHardwareAccelerated && outDim == 24)
            {
                Vector128<float> a0 = Vector128.LoadUnsafe(ref aRef, 0); Vector128<float> a1 = Vector128.LoadUnsafe(ref aRef, 4);
                Vector128<float> a2 = Vector128.LoadUnsafe(ref aRef, 8); Vector128<float> a3 = Vector128.LoadUnsafe(ref aRef, 12);
                Vector128<float> a4 = Vector128.LoadUnsafe(ref aRef, 16); Vector128<float> a5 = Vector128.LoadUnsafe(ref aRef, 20);
                for (int j = 0; j < len; j++)
                {
                    Vector128<float> vx = Vector128.Create(Unsafe.Add(ref xBase, j));
                    ref float wRow = ref Unsafe.Add(ref wBase, j * 24);
                    a0 += vx * Vector128.LoadUnsafe(ref wRow, 0); a1 += vx * Vector128.LoadUnsafe(ref wRow, 4);
                    a2 += vx * Vector128.LoadUnsafe(ref wRow, 8); a3 += vx * Vector128.LoadUnsafe(ref wRow, 12);
                    a4 += vx * Vector128.LoadUnsafe(ref wRow, 16); a5 += vx * Vector128.LoadUnsafe(ref wRow, 20);
                }
                a0.StoreUnsafe(ref aRef, 0); a1.StoreUnsafe(ref aRef, 4); a2.StoreUnsafe(ref aRef, 8);
                a3.StoreUnsafe(ref aRef, 12); a4.StoreUnsafe(ref aRef, 16); a5.StoreUnsafe(ref aRef, 20);
                return;
            }
            if (Vector256.IsHardwareAccelerated && outDim == 32)
            {
                Vector256<float> a0 = Vector256.LoadUnsafe(ref aRef, 0); Vector256<float> a1 = Vector256.LoadUnsafe(ref aRef, 8);
                Vector256<float> a2 = Vector256.LoadUnsafe(ref aRef, 16); Vector256<float> a3 = Vector256.LoadUnsafe(ref aRef, 24);
                for (int j = 0; j < len; j++)
                {
                    Vector256<float> vx = Vector256.Create(Unsafe.Add(ref xBase, j));
                    ref float wRow = ref Unsafe.Add(ref wBase, j * 32);
                    a0 += vx * Vector256.LoadUnsafe(ref wRow, 0); a1 += vx * Vector256.LoadUnsafe(ref wRow, 8);
                    a2 += vx * Vector256.LoadUnsafe(ref wRow, 16); a3 += vx * Vector256.LoadUnsafe(ref wRow, 24);
                }
                a0.StoreUnsafe(ref aRef, 0); a1.StoreUnsafe(ref aRef, 8); a2.StoreUnsafe(ref aRef, 16); a3.StoreUnsafe(ref aRef, 24);
                return;
            }
            if (Vector128.IsHardwareAccelerated && (outDim & 3) == 0)
            {
                for (int j = 0; j < len; j++)
                {
                    Vector128<float> vx = Vector128.Create(Unsafe.Add(ref xBase, j));
                    ref float wRow = ref Unsafe.Add(ref wBase, j * outDim);
                    for (int o = 0; o < outDim; o += 4)
                    {
                        Vector128<float> cur = Vector128.LoadUnsafe(ref aRef, (nuint)o);
                        (cur + vx * Vector128.LoadUnsafe(ref wRow, (nuint)o)).StoreUnsafe(ref aRef, (nuint)o);
                    }
                }
                return;
            }
            for (int j = 0; j < len; j++)
            {
                float xj = Unsafe.Add(ref xBase, j);
                ref float wRow = ref Unsafe.Add(ref wBase, j * outDim);
                for (int o = 0; o < outDim; o++) Unsafe.Add(ref aRef, o) += xj * Unsafe.Add(ref wRow, o);
            }
        }

        private static void ProjectInputs(LinearWeights lin, float[] nodeFeats, int fn, float[] globFeats, int fg, int N, float[] dst)
        {
            int outDim = lin.OutDim;
            Span<float> gContrib = stackalloc float[outDim];
            lin.B.AsSpan().CopyTo(gContrib);
            AccumulateMatVec(lin.WT, fn, globFeats, 0, fg, gContrib, outDim);

            Span<float> acc = stackalloc float[outDim];
            for (int i = 0; i < N; i++)
            {
                gContrib.CopyTo(acc);
                AccumulateMatVec(lin.WT, 0, nodeFeats, i * fn, fn, acc, outDim);
                int yOff = i * outDim;
                for (int o = 0; o < outDim; o++) dst[yOff + o] = acc[o] > 0f ? acc[o] : 0f;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void LinearAddResidualRelu(LinearWeights lin, float[] x, float[] resSrc, int resOff, float[] dst, int dstOff)
        {
            int outDim = lin.OutDim;
            Span<float> acc = stackalloc float[outDim];
            lin.B.AsSpan().CopyTo(acc);
            AccumulateMatVec(lin.WT, 0, x, 0, lin.InDim, acc, outDim);
            for (int o = 0; o < outDim; o++)
            {
                float v = acc[o];
                dst[dstOff + o] = resSrc[resOff + o] + (v > 0f ? v : 0f);
            }
        }

        private static void EvalEdgeHead(LinearWeights h0, LinearWeights h1, float[] hArr, List<Edge> edges, int N, int H, float[] outScores, float[] outSides)
        {
            float[] wt0 = h0.WT, B0 = h0.B, W1 = h1.W;
            float b1 = h1.B[0];
            bool hasSide = h1.OutDim >= 2;
            float b2 = hasSide ? h1.B[1] : 0f;
            Span<float> acc = stackalloc float[H];
            for (int i = 0; i < N; i++)
            {
                if (edges[i].State != EdgeState.Empty) { outScores[i] = -1e9f; outSides[i] = 0f; continue; }
                B0.AsSpan().CopyTo(acc);
                AccumulateMatVec(wt0, 0, hArr, i * H, H, acc, H);
                float score = b1;
                if (hasSide)
                {
                    float side = b2;
                    for (int o = 0; o < H; o++)
                    {
                        float a = acc[o];
                        if (a > 0f) { score += W1[o] * a; side += W1[H + o] * a; }
                    }
                    outSides[i] = side;
                }
                else
                {
                    for (int o = 0; o < H; o++) if (acc[o] > 0f) score += W1[o] * acc[o];
                }
                outScores[i] = score;
            }
        }

        private static void EvalCellHead(LinearWeights h0, LinearWeights h1, float[] hArr, List<Cell> cells, int N, int H, float[] outScores, float[] outSides)
        {
            float[] wt0 = h0.WT, B0 = h0.B, W1 = h1.W;
            float b1 = h1.B[0];
            bool hasSide = h1.OutDim >= 2;
            float b2 = hasSide ? h1.B[1] : 0f;
            Span<float> acc = stackalloc float[H];
            for (int i = 0; i < N; i++)
            {
                if (Math.Abs(cells[i].Color) == 1) { outScores[i] = -1e9f; outSides[i] = 0f; continue; }
                B0.AsSpan().CopyTo(acc);
                AccumulateMatVec(wt0, 0, hArr, i * H, H, acc, H);
                float score = b1;
                if (hasSide)
                {
                    float side = b2;
                    for (int o = 0; o < H; o++)
                    {
                        float a = acc[o];
                        if (a > 0f) { score += W1[o] * a; side += W1[H + o] * a; }
                    }
                    outSides[i] = side;
                }
                else
                {
                    for (int o = 0; o < H; o++) if (acc[o] > 0f) score += W1[o] * acc[o];
                }
                outScores[i] = score;
            }
        }

        private void UpdateDirtyAndClearCooldown(Mesh mesh)
        {
            MeshTopologyCache topo = _topo!;
            int E = topo.NumEdges, C = topo.NumCells;
            List<Edge> edges = mesh.Edges;

            if (!_hasSnapshotState)
            {
                Array.Clear(_edgeFailed, 0, E);
                Array.Clear(_cellFailed, 0, C);
            }
            else
            {
                bool anyRollback = false;
                Array.Clear(_edgeDirty, 0, E);
                for (int e = 0; e < E; e++)
                {
                    byte st = (byte)edges[e].State;
                    if (_lastEdgeState[e] != 0 && st == 0) { anyRollback = true; break; }
                    if (st != _lastEdgeState[e]) _edgeDirty[e] = true;
                }
                if (anyRollback)
                {
                    Array.Clear(_edgeFailed, 0, E);
                    Array.Clear(_cellFailed, 0, C);
                }
                else
                {
                    Array.Copy(_edgeDirty, _edgeDirty2, E);
                    for (int e = 0; e < E; e++)
                    {
                        if (!_edgeDirty[e]) continue;
                        for (int k = topo.EdgeNbrOff[e]; k < topo.EdgeNbrOff[e + 1]; k++)
                            _edgeDirty2[topo.EdgeNbrIdx[k]] = true;
                    }
                    for (int e = 0; e < E; e++)
                    {
                        if (_edgeDirty2[e])
                        {
                            _edgeFailed[e] = false;
                            int c0 = topo.EdgeC0[e], c1 = topo.EdgeC1[e];
                            if (c0 >= 0) _cellFailed[c0] = false;
                            if (c1 >= 0) _cellFailed[c1] = false;
                        }
                    }
                }
            }
            for (int e = 0; e < E; e++) _lastEdgeState[e] = (byte)edges[e].State;
            _hasSnapshotState = true;
        }

        public void OrderIterativeTrials(Mesh mesh, List<IAction> trials, int startIndex, bool isTopLevel, int iterativeRecDepth)
        {
            if (!isTopLevel) return;
            int n = trials.Count - startIndex;
            if (n <= 1) return;

            int minActive = (_autoMode && mesh.IterativeSolverDepth <= 4) ? 24 : _minActiveTrialsToScore;
            int activeCount = 0;
            for (int k = 0; k < n; k++)
                if (!mesh.IsPointlessTrialAction(trials[startIndex + k])) activeCount++;
            if (activeCount < minActive) return;

            EnsureBuffers(mesh);
            int curEmpty = mesh.CountEmptyEdges();
            int cfgHash = ComputeConfigHash(mesh);

            // TwistNTurn bounded-depth coloring trials are very fast (~15-25 us), so Continue + score reuse beats RestartPass.
            // RestartPass + cooldown shines on no-coloring bounded depth (!UseColoring && depth <= 4) and RecDepth >= 2.
            bool useRestartHere = _autoMode
                ? (iterativeRecDepth >= 2 || (!mesh.UseColoring && mesh.IterativeSolverDepth <= 4))
                : (_postAction == PostProgressAction.RestartPass);
            bool useCooldownHere = _autoMode ? useRestartHere : _useSmartCooldown;
            int effRescoreDelta = _autoMode
                ? (iterativeRecDepth >= 2 ? 15 : ((!mesh.UseColoring && mesh.IterativeSolverDepth <= 4) ? 20 : (mesh.IterativeSolverDepth <= 4 ? 12 : 0)))
                : _rescoreEdgeDelta;

            if (iterativeRecDepth != _lastRecDepth || cfgHash != _lastConfigHash)
            {
                _lastRecDepth = iterativeRecDepth;
                _lastConfigHash = cfgHash;
                _hasCachedScores = false;
                _hasSnapshotState = false;
            }

            bool needFullNn = true;
            if (effRescoreDelta > 0 && _hasCachedScores)
            {
                if (curEmpty <= _lastFullNnEmptyEdges && (_lastFullNnEmptyEdges - curEmpty) < effRescoreDelta)
                    needFullNn = false;
            }

            if (needFullNn)
            {
                ScoreAllNodes(mesh, iterativeRecDepth);
                _lastFullNnEmptyEdges = curEmpty;
                _hasCachedScores = true;
            }
            if (useCooldownHere) UpdateDirtyAndClearCooldown(mesh);

            if (_sortScores.Length < n)
            {
                _sortScores = new float[n * 2];
                _sortIndices = new int[n * 2];
                _sortScratch = new IAction[n * 2];
            }

            bool applyFrontierBump = !needFullNn && useCooldownHere;
            bool flipSides = _useSidePrediction && _m.HasSideHead;
            for (int k = 0; k < n; k++)
            {
                IAction action = trials[startIndex + k];
                _sortIndices[k] = k;
                if (mesh.IsPointlessTrialAction(action))
                {
                    _sortScratch[k] = action;
                    _sortScores[k] = -1e9f;
                }
                else if (action is SetAction sa)
                {
                    if (flipSides)
                    {
                        EdgeState desired = _eSides[sa.EdgeIndex] > 0f ? EdgeState.Excluded : EdgeState.Filled;
                        if (sa.EdgeState != desired)
                            action = new SetAction(mesh, sa.EdgeIndex, desired);
                    }
                    _sortScratch[k] = action;
                    float s = _eScores[sa.EdgeIndex];
                    if (applyFrontierBump && _edgeDirty2[sa.EdgeIndex]) s += 1.5f;
                    if (useCooldownHere && _edgeFailed[sa.EdgeIndex]) s -= 100.0f;
                    _sortScores[k] = s;
                }
                else if (action is CellColorJoinAction ccj && ccj.Cell2 == -1)
                {
                    if (flipSides)
                    {
                        bool desiredSame = _cSides[ccj.Cell1] <= 0f;
                        if (ccj.Same != desiredSame)
                            action = new CellColorJoinAction(mesh, ccj.Cell1, -1, desiredSame);
                    }
                    _sortScratch[k] = action;
                    float s = _cScores[ccj.Cell1];
                    if (useCooldownHere && _cellFailed[ccj.Cell1]) s -= 100.0f;
                    _sortScores[k] = s;
                }
                else
                {
                    _sortScratch[k] = action;
                    _sortScores[k] = -1e6f;
                }
            }

            Array.Sort(_sortScores, _sortIndices, 0, n);
            for (int k = 0; k < n; k++)
                trials[startIndex + k] = _sortScratch[_sortIndices[n - 1 - k]];
        }

        public PostProgressAction OnTrialProgress(Mesh mesh, List<IAction> trials, int currentIndex, int addedActions, bool isTopLevel, int iterativeRecDepth)
        {
            if (!isTopLevel) return PostProgressAction.Continue;

            bool useRestartHere = _autoMode
                ? (iterativeRecDepth >= 2 || (!mesh.UseColoring && mesh.IterativeSolverDepth <= 4))
                : (_postAction == PostProgressAction.RestartPass);
            bool useCooldownHere = _autoMode ? useRestartHere : _useSmartCooldown;

            if (useCooldownHere && _topo != null)
            {
                for (int i = 0; i <= currentIndex; i++)
                {
                    IAction a = trials[i];
                    if (a is SetAction sa && sa.EdgeIndex < _edgeFailed.Length) _edgeFailed[sa.EdgeIndex] = true;
                    else if (a is CellColorJoinAction ccj && ccj.Cell2 == -1 && ccj.Cell1 < _cellFailed.Length) _cellFailed[ccj.Cell1] = true;
                }
            }
            if (_autoMode)
                return useRestartHere ? PostProgressAction.RestartPass : PostProgressAction.Continue;

            return addedActions >= _reorderThreshold ? _postAction : PostProgressAction.Continue;
        }

        public int SelectRecursiveTrial(Mesh mesh, List<IAction> trials, int startIndex) => -1;
    }
}
