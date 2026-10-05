using System;
using System.Collections.Generic;
using System.Text;

namespace TwistNTurn
{
    public interface IAction : IEquatable<IAction>
    {
        string Name
        {
            get;
        }
        /// <summary>
        /// Usually true after perform is called if it returns true, but not gauranteed.
        /// </summary>
        bool Successful
        {
            get;
        }
        bool Perform();
        void Unperform();
    }

    public class UndoTree
    {
        public UndoTree()
        {
            root = new UndoNode(null, null, ref curIndex);
            current = root;
        }
        class UndoNode
        {
            public UndoNode(IAction action, UndoNode parent, ref int curIndex)
            {
                Index = curIndex;
                curIndex++;
                Action = action;
                Children = new List<UndoNode>();
                LastDoneChild = null;
                Parent = parent;
                if (parent != null)
                {
                    parent.Children.Add(this);
                    parent.LastDoneChild = this;
                }
                CurrentlyUndone = false;
            }
            public int Index;
            public IAction Action;
            public bool CurrentlyUndone;
            public UndoNode Parent;
            public List<UndoNode> Children;
            public UndoNode LastDoneChild;
        }

        int curIndex = 0;

        UndoNode root;
        UndoNode current;

        public bool Do(IAction action)
        {
            if (action == null)
                return false;
            foreach (UndoNode child in current.Children)
            {
                if (child.Action.Equals(action))
                {
                    current.LastDoneChild = child;
                    return Redo();
                }
            }
            bool res = action.Perform();
            if (!res)
                return false;
            UndoNode node = new UndoNode(action, current, ref curIndex);
            current = node;
            if (markNext)
            {
                marked = current;
                markNext = false;
            }
            OnCanUndoRedoMaybeChanged(EventArgs.Empty);
            return true;
        }

        public event EventHandler CanUndoRedoMaybeChanged;

        private void OnCanUndoRedoMaybeChanged(EventArgs e)
        {
            if (CanUndoRedoMaybeChanged != null)
                CanUndoRedoMaybeChanged(this, e);
        }

        public bool CanUndo
        {
            get
            {
                return current.Action != null;
            }
        }

        public bool Undo()
        {
            if (current.Action == null)
                return false;
            current.Action.Unperform();
            current.CurrentlyUndone = true;
            current = current.Parent;
            OnCanUndoRedoMaybeChanged(EventArgs.Empty);
            return true;
        }

        public bool CanRedo
        {
            get
            {
                return current.LastDoneChild != null;
            }
        }

        public bool Redo()
        {
            UndoNode next = current.LastDoneChild;
            if (next == null)
                return false;
            bool res = next.Action.Perform();
            if (!res)
                return false;
            next.CurrentlyUndone = false;
            current = next;
            if (markNext)
            {
                marked = current;
                markNext = false;
            }
            OnCanUndoRedoMaybeChanged(EventArgs.Empty);
            return true;
        }

        public bool Redo(IAction action)
        {
            if (action == null)
                return false;
            foreach (UndoNode child in current.Children)
            {
                if (child.Action.Equals(action))
                {
                    current.LastDoneChild = child;
                    return Redo();
                }
            }
            return false;
        }

        public List<IAction> CurrentRedoOptions
        {
            get
            {
                List<IAction> res = new List<IAction>();
                foreach (UndoNode child in current.Children)
                    res.Add(child.Action);
                return res;
            }
        }

        internal void Mark()
        {
            marked = current;
        }
        UndoNode marked = null;

        internal void RevertToMark()
        {
            while (current != marked && CanUndo)
                Undo();
        }

        internal object ClearMark()
        {
            UndoNode oldMarked = marked;
            marked = null;
            return oldMarked;
        }

        bool markNext = false;
        internal object MarkNext()
        {
            markNext = true;
            return marked;
        }

        internal void SetMarkedDirect(object prevMarkPos)
        {
            marked = (UndoNode)prevMarkPos;
        }
    }
}
