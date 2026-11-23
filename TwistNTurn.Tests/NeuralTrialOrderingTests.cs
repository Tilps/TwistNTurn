using System;
using TwistNTurn;
using TwistNTurn.TrialOrdering;

namespace TwistNTurn.Tests;

[TestClass]
[DoNotParallelize]
public sealed class NeuralTrialOrderingTests
{
    [TestMethod]
    public void DefaultEmbeddedModel_LoadsWithExpectedMetadata()
    {
        IncidenceGnnModel model = IncidenceGnnModel.Default;
        Assert.AreEqual(2, model.Version);
        Assert.AreEqual(GraphFeatureExtractor.EdgeFeatDim, model.Fe);
        Assert.AreEqual(GraphFeatureExtractor.CellFeatDim, model.Fc);
        Assert.AreEqual(GraphFeatureExtractor.VertFeatDim, model.Fv);
        Assert.AreEqual(GraphFeatureExtractor.GlobalFeatDim, model.Fg);
        Assert.AreEqual(16, model.Hidden);
        Assert.AreEqual(1, model.NumLayers);
        Assert.IsTrue(model.UseColorPool);
        Assert.IsTrue(model.UseMaxPool);
        Assert.IsTrue(model.HasSideHead);
    }

    [TestMethod]
    public void SidePredictionHead_ProducesFiniteSideLogitsAndFlipsPredictedTrials()
    {
        Mesh puzzle = PuzzleHelper.MakeMesh(10, 10, MeshType.Square, 4);
        puzzle.Generate(new Random(101));
        puzzle.Clear();
        puzzle.PerformStart(new System.Collections.Generic.List<IAction>());

        var policy = IncidenceGnnPolicy.CreateDefault();
        policy.ScoreAllNodes(puzzle, iterativeRecDepth: 1);

        int posCount = 0;
        int negCount = 0;
        var trials = new System.Collections.Generic.List<IAction>();
        for (int e = 0; e < puzzle.Edges.Count; e++)
        {
            float logit = policy.EdgeSides[e];
            Assert.IsTrue(float.IsFinite(logit));
            if (puzzle.Edges[e].State == EdgeState.Empty)
            {
                if (logit > 0f) posCount++;
                else if (logit < 0f) negCount++;
                trials.Add(new SetAction(puzzle, e, EdgeState.Filled));
            }
        }

        Assert.IsGreaterThan(0, posCount, "Expected some edges predicted to contradict on Excluded first.");
        Assert.IsGreaterThan(0, negCount, "Expected some edges predicted to contradict on Filled first.");

        policy.OrderIterativeTrials(puzzle, trials, 0, isTopLevel: true, iterativeRecDepth: 1);
        int excludedFirst = 0;
        int filledFirst = 0;
        foreach (IAction action in trials)
        {
            if (action is SetAction sa)
            {
                if (sa.EdgeState == EdgeState.Excluded) excludedFirst++;
                else if (sa.EdgeState == EdgeState.Filled) filledFirst++;
            }
        }

        Assert.IsGreaterThan(0, excludedFirst, "Expected side prediction to flip some trials to Excluded first.");
        Assert.IsGreaterThan(0, filledFirst, "Expected side prediction to keep some trials as Filled first.");
    }

    [TestMethod]
    public void MakeMesh_AttachesPolicyForDifficulties2To5()
    {
        bool prev = PuzzleHelper.UseNeuralTrialOrdering;
        try
        {
            PuzzleHelper.UseNeuralTrialOrdering = true;
            for (int d = 0; d <= 5; d++)
            {
                Mesh mesh = PuzzleHelper.MakeMesh(8, 8, MeshType.Square, d);
                bool expected = d >= 2;
                Assert.AreEqual(expected, mesh.TrialOrderingPolicy != null, $"Difficulty {d}");
            }

            PuzzleHelper.UseNeuralTrialOrdering = false;
            Mesh disabled = PuzzleHelper.MakeMesh(8, 8, MeshType.Square, 4);
            Assert.IsNull(disabled.TrialOrderingPolicy);
        }
        finally
        {
            PuzzleHelper.UseNeuralTrialOrdering = prev;
        }
    }

    [TestMethod]
    [DataRow(MeshType.Square, 10, 10, 4, 101)]
    [DataRow(MeshType.SquareSymmetrical, 10, 10, 4, 202)]
    [DataRow(MeshType.Hexagonal2, 5, 5, 4, 303)]
    [DataRow(MeshType.Hexagonal3, 3, 3, 4, 404)]
    [DataRow(MeshType.Hexagonal4, 3, 3, 4, 505)]
    public void NeuralOrdering_MatchesSequentialSolveAndAnalysis(MeshType type, int w, int h, int diff, int seed)
    {
        bool prev = PuzzleHelper.UseNeuralTrialOrdering;
        try
        {
            PuzzleHelper.UseNeuralTrialOrdering = false;
            Mesh puzzle = PuzzleHelper.MakeMesh(w, h, type, diff);
            puzzle.Generate(new Random(seed));
            puzzle.Clear();

            PuzzleAnalysisResult seqAnalysis = PuzzleHelper.AnalyzePuzzle(puzzle);
            Mesh seqSolver = PuzzleHelper.CreateFullSolverMesh(puzzle);
            seqSolver.Clear();
            SolveState seqState = seqSolver.TrySolve(noRollback: true);

            PuzzleHelper.UseNeuralTrialOrdering = true;
            PuzzleAnalysisResult nnAnalysis = PuzzleHelper.AnalyzePuzzle(puzzle);
            Mesh nnSolver = PuzzleHelper.CreateFullSolverMesh(puzzle);
            nnSolver.Clear();
            SolveState nnState = nnSolver.TrySolve(noRollback: true);

            Assert.AreEqual(SolveState.Solved, seqState);
            Assert.AreEqual(SolveState.Solved, nnState);
            Assert.IsTrue(PuzzleHelper.CheckIsSolved(nnSolver.SolutionFound));
            Assert.AreEqual(seqAnalysis.SolveStatus, nnAnalysis.SolveStatus);
            bool tierMatches = seqAnalysis.TechniqueTier == nnAnalysis.TechniqueTier
                || (seqAnalysis.TechniqueTier.StartsWith("Recursive Search", StringComparison.Ordinal)
                    && nnAnalysis.TechniqueTier.StartsWith("Recursive Search", StringComparison.Ordinal));
            Assert.IsTrue(tierMatches, $"Tier mismatch: {seqAnalysis.TechniqueTier} vs {nnAnalysis.TechniqueTier}");
            Assert.AreEqual(seqAnalysis.BranchDepth, nnAnalysis.BranchDepth);
        }
        finally
        {
            PuzzleHelper.UseNeuralTrialOrdering = prev;
        }
    }
}
