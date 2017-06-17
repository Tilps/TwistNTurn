using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bridge.Html5;
using TwistNTurn;

namespace TwistNTurnBridge
{
    class IntersDisplay
    {

        public IntersDisplay(CanvasRenderingContext2D context, double x, double y)
        {
            this.context = context;
            this.x = x;
            this.y = y;
        }

        private CanvasRenderingContext2D context;
        private double x;
        private double y;
        private bool dirty = true;



        public IntersType IntersType
        {
            get { return intersType; }
            set
            {
                if (value != intersType)
                {
                    intersType = value;
                    dirty = true;
                }
            }
        }

        private IntersType intersType;

        public void Render()
        {
            if (!dirty) return;
            context.GlobalCompositeOperation = CanvasTypes.CanvasCompositeOperationType.DestinationOut;
            context.StrokeStyle = "black";
            context.FillStyle = "black";
            context.BeginPath();
            context.Ellipse(x, y, 10, 10, 0, 0, 360);
            context.Fill();

            context.GlobalCompositeOperation = CanvasTypes.CanvasCompositeOperationType.SourceOver;
            if (IntersType == IntersType.Unknown)
            {
                context.BeginPath();
                context.Ellipse(x, y, 2, 2, 0, 0, 360);
                context.Fill();
            }
            else
            {
                context.BeginPath();
                context.Ellipse(x, y, 8, 8, 0, 0, 360);
                if (IntersType == IntersType.Twist)
                {
                    context.FillStyle = "white";
                }
                context.Fill();
                context.Stroke();

            }
        }

    }
}
