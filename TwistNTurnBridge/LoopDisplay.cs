using System;
using System.Collections.Generic;
using System.Linq;
using Bridge.Html5;
using TwistNTurn;


namespace TwistNTurnBridge
{
    public class LoopDisplay
    {
        public LoopDisplay(HTMLDivElement displayHost)
        {
            this.displayHost = displayHost;
            this.Mesh = new Mesh(10, 10, MeshType.Square);
            displayHost.OnMouseDown = LoopDisplay_MouseDown;
            displayHost.OnMouseUp = LoopDisplay_MouseUp;
            displayHost.OnMouseLeave = LoopDisplay_MouseLeave;
            displayHost.OnMouseMove = LoopDisplay_MouseMove;
            displayHost.OnContextMenu = LoopDisplay_ContextMenu;
            cellColorCanvas = Document.CreateElement<HTMLCanvasElement>("canvas");
            cellColorCanvas.Style.Position = Position.Absolute;
            cellColorCanvas.Style.ZIndex = "1";
            cellNumberCanvas = Document.CreateElement<HTMLCanvasElement>("canvas");
            cellNumberCanvas.Style.Position = Position.Absolute;
            cellNumberCanvas.Style.ZIndex = "2";
            edgeCanvas = Document.CreateElement<HTMLCanvasElement>("canvas");
            edgeCanvas.Style.Position = Position.Absolute;
            edgeCanvas.Style.ZIndex = "3";
            intersectionCanvas = Document.CreateElement<HTMLCanvasElement>("canvas");
            intersectionCanvas.Style.Position = Position.Absolute;
            intersectionCanvas.Style.ZIndex = "4";
            displayHost.AppendChild(cellColorCanvas);
            displayHost.AppendChild(cellNumberCanvas);
            displayHost.AppendChild(edgeCanvas);
            displayHost.AppendChild(intersectionCanvas);
        }

        private HTMLDivElement displayHost;
        private HTMLCanvasElement cellColorCanvas;
        private HTMLCanvasElement cellNumberCanvas;
        private HTMLCanvasElement edgeCanvas;
        private HTMLCanvasElement intersectionCanvas;

        void LoopDisplay_ContextMenu(Event e)
        {
            clicking = false;
            e.PreventDefault();
            e.StopPropagation();
        }

        void LoopDisplay_MouseLeave(MouseEvent e)
        {
            clicking = false;
        }

        bool clicking = false;

        void LoopDisplay_MouseUp(MouseEvent e)
        {
            e.PreventDefault();
            e.StopPropagation();

            clicking = false;
        }

        void LoopDisplay_MouseMove(MouseEvent e)
        {
            if (clicking)
            {
                int closestEdge;
                int closestCell;
                var bounds = displayHost.GetBoundingClientRect();
                FindClosest(e.ClientX - bounds.Left , e.ClientY - bounds.Top, out closestEdge, out closestCell, 8.0f / (float)scaleFactor);
                if (closestEdge != -1 || closestCell != -1)
                {
                    if (lastClosestCell != closestCell || lastClosestEdge != closestEdge)
                    {
                        lastClosestEdge = closestEdge;
                        lastClosestCell = closestCell;
                        PerformAction(null, closestEdge, closestCell);
                    }
                }
            }
        }

        void LoopDisplay_MouseDown(MouseEvent e)
        {
            e.PreventDefault();
            e.StopPropagation();
            clicking = true;
            OnMouseButtonDown(e);
        }

        HashSet<int> markedEdges = new HashSet<int>();

        internal void OnKeyDown(KeyboardEvent e)
        {

            if (e.CtrlKey && e.KeyCode == KeyboardEvent.DOM_VK_Z)
            {
                if (this.undoTree.CanUndo)
                {
                    this.undoTree.Undo();
                    this.UpdateChildControls();
                }
                e.PreventDefault();
                e.StopPropagation();
            }
            else if (e.CtrlKey && e.KeyCode == KeyboardEvent.DOM_VK_Y)
            {
                if (this.undoTree.CanRedo)
                {
                    this.undoTree.Redo();
                    this.UpdateChildControls();
                }
                e.PreventDefault();
                e.StopPropagation();
            }
            else if (e.CtrlKey && e.KeyCode == KeyboardEvent.DOM_VK_F)
            {
                FixPosition();
                e.PreventDefault();
                e.StopPropagation();
            }
            else if (e.CtrlKey && (e.KeyCode == KeyboardEvent.DOM_VK_R || e.KeyCode == KeyboardEvent.DOM_VK_M))
            {
                ResetToFixed();
                e.PreventDefault();
                e.StopPropagation();
            }
            else if (e.CtrlKey && e.KeyCode == KeyboardEvent.DOM_VK_U)
            {
                Unfix();
                e.PreventDefault();
                e.StopPropagation();
            }
        }

        internal void Unfix()
        {
            this.undoTree.Do(new FixAction(this, false));
            this.UpdateChildControls();
        }

        private void ResetToFixed()
        {
            this.undoTree.RevertToMark();
            this.UpdateChildControls();
        }

        private void FixPosition()
        {
            this.undoTree.Do(new FixAction(this, true));
            this.UpdateChildControls();
        }


        public Mesh Mesh
        {
            get
            {
                return mesh;
            }
            set
            {
                this.mesh = value;
                markedEdges.Clear();
                undoTree = new UndoTree();
                OnMeshChanged(mesh);
            }
        }

        private Mesh mesh;


        private double scaleFactor = 20;
        private double xOffset = 10;
        private double yOffset = 10;

        internal void OnScaleChanged()
        {
            OnMeshChanged(this.Mesh);
        }

        private void UpdateChildControls()
        {
            for (int i = 0; i < this.Mesh.Edges.Count; i++)
            {
                Edge edge = this.Mesh.Edges[i];
                EdgeDisplay edgeDisplay = edgeDisplays[i];
                edgeDisplay.EdgeState = edge.State;
                //edgeDisplay.EdgeColor = edge.Color;
                edgeDisplay.Marked = markedEdges.Contains(i);
                edgeDisplay.Render();
            }
        }
        List<EdgeDisplay> edgeDisplays = new List<EdgeDisplay>();
        List<IntersDisplay> intersDisplays = new List<IntersDisplay>();

        private void OnMeshChanged(Mesh newMesh)
        {
            if (this.displayHost.ClientHeight <= 0 || this.displayHost.ClientWidth <= 0)
                return;
            UpdateScaleFactor(newMesh);
            edgeDisplays.Clear();
            intersDisplays.Clear();
            if (newMesh == null)
                return;
            var drawContext = edgeCanvas.GetContext(CanvasTypes.CanvasContext2DType.CanvasRenderingContext2D);
            drawContext.ClearRect(0, 0, edgeCanvas.Width, edgeCanvas.Height);
            for (int i = 0; i < newMesh.Edges.Count; i++)
            {
                Edge edge = newMesh.Edges[i];
                float x1, x2, y1, y2;
                newMesh.GetEdgeExtent(edge, out x1, out y1, out x2, out y2);
                EdgeDisplay edgeDisplay = new EdgeDisplay(drawContext, x1*scaleFactor+xOffset, x2*scaleFactor+xOffset, y1*scaleFactor+yOffset, y2*scaleFactor+yOffset);
                edgeDisplay.EdgeState = edge.State;
                //edgeDisplay.EdgeColor = edge.Color;
                edgeDisplay.Marked = markedEdges.Contains(i);
                edgeDisplay.Render();
                edgeDisplays.Add(edgeDisplay);
            }

            drawContext = intersectionCanvas.GetContext(CanvasTypes.CanvasContext2DType.CanvasRenderingContext2D);
            drawContext.ClearRect(0, 0, intersectionCanvas.Width, intersectionCanvas.Height);
            foreach (var intersection in mesh.Intersections)
            {
                IntersDisplay intersDisplay = new IntersDisplay(drawContext, intersection.X * scaleFactor + xOffset, intersection.Y * scaleFactor + yOffset);
                intersDisplay.IntersType = intersection.Type;
                intersDisplay.Render();
                intersDisplays.Add(intersDisplay);
            }
        }

        private double emRatio = -1;

        private void UpdateScaleFactor(Mesh newMesh)
        {
            if (newMesh == null)
                return;
            cellColorCanvas.Width = displayHost.ClientWidth;
            cellColorCanvas.Height = displayHost.ClientHeight;
            cellNumberCanvas.Width = displayHost.ClientWidth;
            cellNumberCanvas.Height = displayHost.ClientHeight;
            edgeCanvas.Width = displayHost.ClientWidth;
            edgeCanvas.Height = displayHost.ClientHeight;
            intersectionCanvas.Width = displayHost.ClientWidth;
            intersectionCanvas.Height = displayHost.ClientHeight;
            double minX = double.MaxValue;
            double minY = double.MaxValue;
            double maxX = double.MinValue;
            double maxY = double.MinValue;
            for (int i = 0; i < newMesh.Intersections.Count; i++)
            {
                Intersection inters = newMesh.Intersections[i];
                if (inters.X < minX)
                    minX = inters.X;
                if (inters.Y < minY)
                    minY = inters.Y;
                if (inters.X > maxX)
                    maxX = inters.X;
                if (inters.Y > maxY)
                    maxY = inters.Y;
            }
            double marginTop = 2;
            double marginBottom = 2;
            double marginLeft = 2;
            double marginRight = 2;
            scaleFactor = (displayHost.ClientHeight - marginTop - marginBottom) / (maxY - minY + 1);
            scaleFactor = Math.Min(scaleFactor, (displayHost.ClientWidth - marginLeft - marginRight) / (maxX - minX + 1));
            scaleFactor = Math.Max(scaleFactor, 25.0);
            scaleFactor = Math.Min(scaleFactor, 50.0);
            yOffset = displayHost.ClientHeight / 2.0 - (maxY + minY) * scaleFactor / 2.0;
            xOffset = displayHost.ClientWidth / 2.0 - (maxX + minX) * scaleFactor / 2.0;
        }

        int lastControl = -1;
        int lastShift = -1;


        private double DistanceToSegmentSq(double realX, double realY, double sx, double sy, double ex, double ey)
        {
            double dx = ex - sx;
            double dy = ey - sy;
            double t = ((realX - sx) * dx + (realY - sy) * dy) / (dx * dx + dy * dy);
            double nearX;
            double nearY;
            if (t < 0)
            {
                nearX = sx;
                nearY = sy;
            }
            else if (t > 1)
            {
                nearX = ex;
                nearY = ey;
            }
            else
            {
                nearX = sx + t * dx;
                nearY = sy + t * dy;
            }
            double distx = nearX - realX;
            double disty = nearY - realY;
            return distx * distx + disty * disty;
        }

        private bool showCellColors = false;
        public bool ShowColors
        {
            get
            {
                return showColors;
            }
            set
            {
                showColors = value;
            }
        }
        private bool showColors = false;
        private bool noToggle = false;
        private int autoMove = 0;
        private bool disallowFalseMove = false;
        private bool useICInAuto = false;
        private bool considerMultipleLoopsInAuto = false;
        private bool useCellColoringInAuto = false;
        private bool useColoringInAuto = false;


        private UndoTree undoTree = new UndoTree();


        internal void OnMouseButtonDown(MouseEvent e)
        {
            if (e.Button != 0)
            {
                e.PreventDefault();
                e.StopPropagation();
            }
            var bounds = displayHost.GetBoundingClientRect();
            OnMouseButtonDown(e.ClientX - bounds.Left, e.ClientY - bounds.Top, e.Button != 0, e.ShiftKey, e.CtrlKey);
        }
        internal void OnMouseButtonDown(double clickX, double clickY, bool right, bool shiftPressed, bool controlPressed)
        {
            if (!ShowColors && shiftPressed)
                right = !right;
            else if (!ShowColors && (shiftPressed || controlPressed))
                return;
            if (ShowColors && !shiftPressed)
                lastShift = -1;
            if (ShowColors && !controlPressed)
                lastControl = -1;
            int closestEdge;
            int closestCell;
            FindClosest(clickX, clickY, out closestEdge, out closestCell, 0.0f);
            if (closestEdge != -1 || closestCell != -1)
            {
                lastClosestEdge = closestEdge;
                lastClosestCell = closestCell;
            }
            PerformAction(right, closestEdge, closestCell);
        }
        int lastClosestEdge = -1;
        int lastClosestCell = -1;

        private void FindClosest(double clickX, double clickY, out int closestEdge, out int closestCell, float minSep)
        {
            double linkLength = scaleFactor;
            double x = clickX - xOffset;
            double y = clickY - yOffset;
            double realX = x / linkLength;
            double realY = y / linkLength;
            closestEdge = -1;
            closestCell = -1;
            double distanceSq = 1;
            double secondDistanceSq = 200*200;
            for (int i = 0; i < this.Mesh.Edges.Count; i++)
            {
                Edge edge = this.Mesh.Edges[i];
                float sx, sy, ex, ey;
                this.Mesh.GetEdgeExtent(edge, out sx, out sy, out ex, out ey);
                double curDistSq = DistanceToSegmentSq(realX, realY, sx, sy, ex, ey);
                if (curDistSq <= distanceSq)
                {
                    closestEdge = i;
                    secondDistanceSq = distanceSq;
                    distanceSq = curDistSq;
                }
                else if (curDistSq <= secondDistanceSq)
                {
                    secondDistanceSq = curDistSq;
                }
            }
            if (showCellColors)
            {
                for (int i = 0; i < this.Mesh.Cells.Count; i++)
                {
                    float cx = 0.0F;
                    float cy = 0.0F;
                    foreach (int inters in Mesh.Cells[i].Intersections)
                    {
                        cx += Mesh.Intersections[inters].X;
                        cy += Mesh.Intersections[inters].Y;
                    }
                    cx /= Mesh.Cells[i].Intersections.Count;
                    cy /= Mesh.Cells[i].Intersections.Count;
                    double distx = cx - realX;
                    double disty = cy - realY;
                    double curDistSq = distx * distx + disty * disty;
                    if (curDistSq <= distanceSq)
                    {
                        closestEdge = -1;
                        closestCell = i;
                        secondDistanceSq = distanceSq;
                        distanceSq = curDistSq;
                    }
                    else if (curDistSq <= secondDistanceSq)
                    {
                        secondDistanceSq = curDistSq;
                    }
                }
            }
            if (distanceSq + 2*minSep*Math.Sqrt(distanceSq) + minSep*minSep > secondDistanceSq)
            {
                closestEdge = -1;
                closestCell = -1;
            }
        }


        EdgeState lastState = EdgeState.Filled;

        private void PerformAction(bool? right, int closestEdge, int closestCell)
        {
            if (closestEdge != -1)
            {
                if (markedEdges.Contains(closestEdge))
                    return;
                if (noToggle && Mesh.Edges[closestEdge].State != EdgeState.Empty)
                    return;
                /*if (shiftPressed || controlPressed)
                {
                    if (shiftPressed)
                    {
                        if (lastShift == -1)
                            lastShift = closestEdge;
                        else
                        {
                            ColorJoinAction colorAction = new ColorJoinAction(Mesh, lastShift, closestEdge, true);
                            if (lastShift != closestEdge)
                                undoTree.Do(colorAction);
                            lastShift = -1;
                        }
                    }
                    else if (controlPressed)
                    {
                        if (lastControl == -1)
                            lastControl = closestEdge;
                        else
                        {
                            ColorJoinAction colorAction = new ColorJoinAction(Mesh, lastControl, closestEdge, false);
                            if (lastControl != closestEdge)
                                undoTree.Do(colorAction);
                            lastControl = -1;
                        }
                    }
                    UpdateChildControls();
                    return;
                }*/
                LoopClickAction action;
                if (right.HasValue)
                {
                    action = new LoopClickAction(Mesh, closestEdge, right.Value, autoMove, disallowFalseMove, useICInAuto, considerMultipleLoopsInAuto, useColoringInAuto, useCellColoringInAuto);
                }
                else
                {
                    if (Mesh.Edges[closestEdge].State == lastState)
                        return;
                    bool pretendRight = false;
                    switch (lastState)
                    {
                        case EdgeState.Filled:
                            if (Mesh.Edges[closestEdge].State == EdgeState.Excluded)
                                pretendRight = true;
                            break;
                        case EdgeState.Excluded:
                            if (Mesh.Edges[closestEdge].State == EdgeState.Empty)
                                pretendRight = true;
                            break;
                        case EdgeState.Empty:
                            if (Mesh.Edges[closestEdge].State == EdgeState.Filled)
                                pretendRight = true;
                            break;
                    }
                    action = new LoopClickAction(Mesh, closestEdge, pretendRight, autoMove, disallowFalseMove, useICInAuto, considerMultipleLoopsInAuto, useColoringInAuto, useCellColoringInAuto);
                }
                if (!undoTree.Do(action))
                {
                    /*
                    redEdge = closestEdge;
                    Thread thread = new Thread(new ThreadStart(ClearRed));
                    thread.IsBackground = true;
                    thread.Start();
                     */
                }
                else if (noToggle)
                {
                    /*
                    if (MovePerformed != null)
                        MovePerformed(this, new MoveEventArgs(closestEdge, e.Button == MouseButtons.Left));
                     * */
                }
                else
                {
                    if (right.HasValue)
                    {
                        lastState = Mesh.Edges[closestEdge].State;
                    }
                    bool satisified = true;
                    bool nonempty = false;
                    for (int i = 0; i < Mesh.Intersections.Count; i++)
                    {
                        if (Mesh.Intersections[i].FilledCount != 2 && Mesh.Intersections[i].FilledCount != 0)
                        {
                            satisified = false;
                            break;
                        }
                        if (Mesh.Intersections[i].Type != IntersType.Unknown)
                        {
                            nonempty = true;
                            if (!Mesh.TypeSatisfied(i))
                            {
                                satisified = false;
                                break;
                            }
                        }
                    }
                    if (satisified && nonempty)
                    {
                        Mesh copy = new Mesh(Mesh);
                        try
                        {
                            copy.Clear();
                            bool failed = false;
                            if (copy.TrySolve() != SolveState.Solved)
                            {
                                copy.SolverMethod = SolverMethod.Recursive;
                                //copy.UseIntersectCellInteractsInSolver = false;
                                copy.UseCellColoringTrials = false;
                                copy.UseCellColoring = true;
                                copy.UseCellPairs = false;
                                copy.UseCellPairsTopLevel = true;
                                copy.UseColoring = true;
                                copy.UseDerivedColoring = true;
                                copy.UseEdgeRestricts = true;
                                copy.UseMerging = true;
                                copy.ConsiderMultipleLoops = true;
                                copy.ColoringCheats = true;
                                if (copy.TrySolve() != SolveState.Solved)
                                {
                                    failed = true;
                                }
                            }
                            if (!failed)
                            {
                                bool done = true;
                                for (int i = 0; i < Mesh.Edges.Count; i++)
                                {
                                    if (copy.SolutionFound.Edges[i].State == EdgeState.Filled)
                                    {
                                        if (Mesh.Edges[i].State != EdgeState.Filled)
                                            done = false;
                                    }
                                    else if (Mesh.Edges[i].State == EdgeState.Filled)
                                        done = false;
                                }
                                if (done)
                                {
                                    FixPosition();
                                    copy.FullClear();
                                }
                            }
                        }
                        catch
                        {
                        }
                        
                    }
                }
                UpdateChildControls();
            }
            else if (closestCell != -1)
            {
                if (noToggle && Mesh.Cells[closestCell].Color != 0)
                    return;
                CellClickAction action = new CellClickAction(Mesh, closestCell, right.Value);
                undoTree.Do(action);
                UpdateChildControls();
            }
        }

        internal List<int> FixPositionInternal(bool fix, out object lastMark)
        {
            List<int> result = markedEdges.ToList();
            markedEdges.Clear();
            if (fix)
            {
                for (int i = 0; i < this.Mesh.Edges.Count; i++)
                {
                    if (this.Mesh.Edges[i].State != EdgeState.Empty)
                        markedEdges.Add(i);
                }
                lastMark = this.undoTree.MarkNext();
            }
            else
            {
                lastMark = this.undoTree.ClearMark();
            }
            return result;
        }

        internal void SetMarkedDirectInternal(List<int> prevMarked)
        {
            markedEdges.Clear();
            markedEdges.UnionWith(prevMarked);
        }

        internal void SetPrevMarkPos(object prevMarkPos)
        {
            this.undoTree.SetMarkedDirect(prevMarkPos);
        }

        internal void Undo()
        {
            if (undoTree.CanUndo)
            {
                undoTree.Undo();
                UpdateChildControls();
            }
        }

        internal void Redo()
        {
            if (undoTree.CanRedo)
            {
                undoTree.Redo();
                UpdateChildControls();
            }
        }


        internal void Fix()
        {
            FixPosition();
        }

        internal void RevertToFix()
        {
            ResetToFixed();
        }
    }

    class FixAction : IAction
    {
        public FixAction(LoopDisplay display, bool fix)
        {
            this.display = display;
            this.fix = fix;
        }
        LoopDisplay display;
        bool fix = true;
        List<int> prevMarked;
        object prevMarkPos;

        public string Name
        {
            get { return fix ? "Fix Position" : "Unfix Position"; }
        }

        public bool Successful
        {
            get { return true; }
        }

        public bool Perform()
        {
            prevMarked = display.FixPositionInternal(fix, out prevMarkPos);
            return true;
        }

        public void Unperform()
        {
            display.SetMarkedDirectInternal(prevMarked);
            display.SetPrevMarkPos(prevMarkPos);
        }

        public bool Equals(IAction other)
        {
            FixAction realOther = other as FixAction;
            if (realOther == null)
                return false;
            return realOther.display == this.display && realOther.fix == this.fix;
        }
    }

    class LoopClickAction : IAction
    {
        public LoopClickAction(Mesh mesh, int edgeIndex, bool buttons, int autoMove, bool disallowFalseMove, bool useICInAuto, bool considerMultipleLoopsInAuto, bool useColoringInAuto, bool useCellColoringInAuto)
        {
            this.mesh = mesh;
            this.edgeIndex = edgeIndex;
            this.buttons = buttons;
            this.autoMove = autoMove;
            this.disallowFalseMove = disallowFalseMove;
            this.useICInAuto = useICInAuto;
            this.considerMultipleLoopsInAuto = considerMultipleLoopsInAuto;
            this.useColoringInAuto = useColoringInAuto;
            this.useCellColoringInAuto = useCellColoringInAuto;
        }

        Mesh mesh;
        int edgeIndex;
        bool buttons;
        int autoMove;
        bool disallowFalseMove;
        bool useICInAuto;
        bool considerMultipleLoopsInAuto;
        bool useColoringInAuto;
        bool useCellColoringInAuto;

        List<IAction> actionsPerformed;

        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;



        private EdgeState Toggle(EdgeState loopLinkState, bool mouseButtons)
        {
            if (!mouseButtons)
            {
                if (loopLinkState == EdgeState.Filled)
                    return EdgeState.Excluded;
                else if (loopLinkState == EdgeState.Excluded)
                    return EdgeState.Empty;
                else
                    return EdgeState.Filled;
            }
            else
            {
                if (loopLinkState == EdgeState.Excluded)
                    return EdgeState.Filled;
                else if (loopLinkState == EdgeState.Filled)
                    return EdgeState.Empty;
                else
                    return EdgeState.Excluded;
            }
        }
        #region IAction Members

        public string Name
        {
            get
            {
                string clickName = string.Empty;
                if (!buttons)
                    clickName = "Left Click";
                else
                    clickName = "Right Click";
                return clickName + " Edge: " + edgeIndex.ToString();
            }
        }

        public bool Perform()
        {
            successful = true;
            // mesh.ConsiderIntersectCellInteractsAsSimple = useICInAuto;
            mesh.ConsiderMultipleLoops = considerMultipleLoopsInAuto;
            mesh.UseColoring = useColoringInAuto;
            mesh.UseCellColoring = useCellColoringInAuto;
            // TODO: add settings for use derived/merge in auto.
            mesh.UseDerivedColoring = false;
            mesh.UseMerging = false;
            mesh.ColoringCheats = false;
            Edge closest = mesh.Edges[edgeIndex];
            EdgeState toggled = Toggle(closest.State, buttons);
            actionsPerformed = new List<IAction>();
            if (closest.State != EdgeState.Empty)
            {
                IAction unsetAction = new UnsetAction(mesh, edgeIndex);
                bool res = unsetAction.Perform();
                if ((!res || !unsetAction.Successful) && disallowFalseMove)
                {
                    if (res && !unsetAction.Successful)
                        actionsPerformed.Add(unsetAction);
                    Unperform();
                    return false;
                }
                else if (!res || !unsetAction.Successful)
                    successful = false;
                actionsPerformed.Add(unsetAction);
            }
            if (toggled != EdgeState.Empty)
            {
                bool res = mesh.Perform(edgeIndex, toggled, actionsPerformed, autoMove);
                if (!res && disallowFalseMove)
                {
                    Unperform();
                    return false;
                }
                else if (!res)
                    successful = false;
            }
            return true;
        }

        public void Unperform()
        {
            if (actionsPerformed.Count > 0)
            {
                mesh.Unperform(actionsPerformed);
            }
        }

        #endregion

        #region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            LoopClickAction realOther = other as LoopClickAction;
            if (realOther == null)
                return false;
            if (this.mesh == realOther.mesh &&
                realOther.buttons == this.buttons &&
                realOther.edgeIndex == this.edgeIndex &&
                this.autoMove == realOther.autoMove &&
                this.considerMultipleLoopsInAuto == realOther.considerMultipleLoopsInAuto &&
                this.useICInAuto == realOther.useICInAuto &&
                this.useColoringInAuto == realOther.useColoringInAuto &&
                this.disallowFalseMove == realOther.disallowFalseMove)
                return true;
            return false;
        }

        #endregion
    }

    class CellClickAction : IAction
    {
        public CellClickAction(Mesh mesh, int cellIndex, bool buttons)
        {
            this.mesh = mesh;
            this.cellIndex = cellIndex;
            this.buttons = buttons;
        }

        Mesh mesh;
        int cellIndex;
        bool buttons;

        List<IAction> actionsPerformed;

        public bool Successful
        {
            get
            {
                return successful;
            }
        }
        private bool successful;



        private int Toggle(int color, bool mouseButtons)
        {
            if (!mouseButtons)
            {
                if (color == 1)
                    return -1;
                else if (color == -1)
                    return 0;
                return 1;
            }
            else
            {
                if (color == 1)
                    return 0;
                else if (color == -1)
                    return 1;
                return -1;
            }
        }
        #region IAction Members

        public string Name
        {
            get
            {
                string clickName = string.Empty;
                if (!buttons)
                    clickName = "Left Click";
                else
                    clickName = "Right Click";
                return clickName + " Cell: " + cellIndex.ToString();
            }
        }

        public bool Perform()
        {
            successful = true;
            Cell closest = mesh.Cells[cellIndex];
            int newColor = Toggle(closest.Color, buttons);
            actionsPerformed = new List<IAction>();
            if (closest.Color != 0)
            {
                IAction unsetAction = new CellColorClearAction(mesh, cellIndex);
                bool res = unsetAction.Perform();
                if ((!res || !unsetAction.Successful))
                {
                    if (res && !unsetAction.Successful)
                        actionsPerformed.Add(unsetAction);
                    Unperform();
                    return false;
                }
                actionsPerformed.Add(unsetAction);
            }
            if (newColor != 0)
            {
                IAction setAction = new CellColorJoinAction(mesh, cellIndex, -1, newColor == 1);
                bool res = setAction.Perform();
                if ((!res || !setAction.Successful))
                {
                    if (res && !setAction.Successful)
                        actionsPerformed.Add(setAction);
                    Unperform();
                    return false;
                }
                actionsPerformed.Add(setAction);
            }
            return true;
        }

        public void Unperform()
        {
            if (actionsPerformed.Count > 0)
            {
                mesh.Unperform(actionsPerformed);
            }
        }

        #endregion

        #region IEquatable<IAction> Members

        public bool Equals(IAction other)
        {
            CellClickAction realOther = other as CellClickAction;
            if (realOther == null)
                return false;
            if (this.mesh == realOther.mesh &&
                realOther.buttons == this.buttons &&
                realOther.cellIndex == this.cellIndex)
                return true;
            return false;
        }

        #endregion
    }
}
