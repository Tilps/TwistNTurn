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
            public UndoNode(IAction? action, UndoNode? parent, ref int curIndex)
            {
                Index = curIndex;
                curIndex++;
                Action = action;
                Children = [];
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
            public IAction? Action;
            public bool CurrentlyUndone;
            public UndoNode? Parent;
            public List<UndoNode> Children;
            public UndoNode? LastDoneChild;
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
                if (child.Action != null && child.Action.Equals(action))
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

        public event EventHandler? CanUndoRedoMaybeChanged;

        private void OnCanUndoRedoMaybeChanged(EventArgs e)
        {
            CanUndoRedoMaybeChanged?.Invoke(this, e);
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
            if (current.Action == null || current.Parent == null)
                return false;
            current.Action.Unperform();
            current.CurrentlyUndone = true;
            current = current.Parent;
            OnCanUndoRedoMaybeChanged(EventArgs.Empty);
            return true;
        }

        public bool UndoAndForget()
        {
            if (current.Action == null || current.Parent == null)
                return false;
            var toRemove = current;
            current.Action.Unperform();
            current = current.Parent;
            current.Children.Remove(toRemove);
            if (current.LastDoneChild == toRemove)
            {
                current.LastDoneChild = current.Children.Count > 0 ? current.Children[^1] : null;
            }
            if (marked == toRemove)
            {
                marked = current;
            }
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
            UndoNode? next = current.LastDoneChild;
            if (next == null || next.Action == null)
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
                if (child.Action != null && child.Action.Equals(action))
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
                List<IAction> res = [with(current.Children.Count)];
                foreach (UndoNode child in current.Children)
                {
                    if (child.Action != null)
                        res.Add(child.Action);
                }
                return res;
            }
        }

        public void Mark()
        {
            marked = current;
        }
        UndoNode? marked = null;

        public bool IsMarked => marked != null;
        public bool CanRevert => marked != null && current != marked;

        public void RevertToMark()
        {
            while (current != marked && CanUndo)
                Undo();
        }

        public object? ClearMark()
        {
            UndoNode? oldMarked = marked;
            marked = null;
            return oldMarked;
        }

        bool markNext = false;
        public object? MarkNext()
        {
            markNext = true;
            return marked;
        }

        public void SetMarkedDirect(object? prevMarkPos)
        {
            marked = prevMarkPos as UndoNode;
        }
    }
}
