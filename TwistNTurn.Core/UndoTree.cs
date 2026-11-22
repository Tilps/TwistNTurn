using System;
using System.Collections.Generic;

namespace TwistNTurn
{
    public interface IAction : IEquatable<IAction>
    {
        string Name
        {
            get;
        }
        /// <summary>
        /// Usually true after perform is called if it returns true, but not guaranteed.
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
            root = new UndoNode(null, null);
            current = root;
        }

        class UndoNode
        {
            public UndoNode(IAction? action, UndoNode? parent)
            {
                Action = action;
                Children = [];
                LastDoneChild = null;
                Parent = parent;
                if (parent != null)
                {
                    parent.Children.Add(this);
                    parent.LastDoneChild = this;
                }
            }
            public IAction? Action;
            public UndoNode? Parent;
            public List<UndoNode> Children;
            public UndoNode? LastDoneChild;
        }

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
                    return RedoChild(child);
                }
            }
            bool res = action.Perform();
            if (!res)
                return false;
            UndoNode node = new UndoNode(action, current);
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

        public bool CanUndo => current.Action != null;

        public IAction? CurrentAction => current.Action;

        public bool Undo()
        {
            if (current.Action == null || current.Parent == null)
                return false;
            current.Action.Unperform();
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
            if (marked != null && IsDescendantOrSelf(marked, toRemove))
            {
                marked = current;
            }
            OnCanUndoRedoMaybeChanged(EventArgs.Empty);
            return true;
        }

        public bool CanRedo => current.LastDoneChild != null;

        public bool Redo()
        {
            UndoNode? next = current.LastDoneChild;
            if (next == null)
                return false;
            return RedoChild(next);
        }

        private bool RedoChild(UndoNode next)
        {
            if (next.Action == null)
                return false;
            bool res = next.Action.Perform();
            if (!res)
                return false;
            current.LastDoneChild = next;
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
                    return RedoChild(child);
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
        public bool CanRevert => marked != null && current != marked && IsDescendantOrSelf(current, marked);

        private static bool IsDescendantOrSelf(UndoNode node, UndoNode potentialAncestor)
        {
            UndoNode? cursor = node;
            while (cursor != null)
            {
                if (cursor == potentialAncestor)
                    return true;
                cursor = cursor.Parent;
            }
            return false;
        }

        public void RevertToMark()
        {
            if (marked == null || !IsDescendantOrSelf(current, marked))
                return;
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

        public object? GetMarkedDirect() => marked;

        public void SetMarkedDirect(object? prevMarkPos)
        {
            marked = prevMarkPos as UndoNode;
        }
    }
}
