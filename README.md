# Twist 'n' Turn

**Twist 'n' Turn** is a logic puzzle game based on **Masyu** (also known as Shiroshinju Kuroshinju / White Pearls and Black Pearls), played on planar grids. The goal is to draw a single continuous, non-intersecting loop along the grid lines that passes through all pearl clues according to their rules.

The project features a solver and generator that work across arbitrary planar mesh topologies (including Square, Square Symmetrical, and Hexagonal grids), running in the browser via a modern Blazor WebAssembly frontend.

## Pearl Rules

- **⚪ White Pearl (Twist):** The loop must pass *straight* through the pearl's intersection, and must turn 90° in at least one of the intersections immediately before or after it.
- **⚫ Black Pearl (Turn):** The loop must make a 90° *turn* at the pearl's intersection, and must travel straight through the intersections immediately before and after it.
- **Grid Intersections:** Empty grid intersections without pearls may be traversed straight, turned at, or not visited at all. Every visited intersection must have exactly 2 loop edges attached (no branches or crossings).

## Grid Topologies Supported

- **Square**
- **Square Symmetrical**
- **Hexagon 3**

## Architecture

The project is built on **.NET 11** and **C# 15**, organized into modular SDK-style projects with zero compiler warnings and nullable reference types enabled:

- **`TwistNTurn.Core`** (`net11.0`): Core puzzle generation and solver engine supporting arbitrary planar graphs and depth-limited deductive inference, URL-safe compact codec/compression, and an undo/redo action tree.
- **`TwistNTurn.Web`** (`net11.0` Blazor WebAssembly): Modern, responsive web frontend running client-side in the browser via WebAssembly (with AOT compilation in Release). Features SVG rendering for Masyu white and black pearls, single-finger drawing/cycling, two-finger pinch-to-zoom and pan, keyboard shortcuts, checkpoints (Fix/Unfix/Revert), URL hash sharing, and PWA offline support.
- **`TwistNTurn.Tests`** (`net11.0` MSTest): Automated test suite covering grid topology creation, pearl validation rules, solver deductions, undo/redo trees, edge state serialization, and puzzle codec round-tripping.

## Building and Running

### Prerequisites
- [.NET 11 SDK](https://dotnet.microsoft.com/)
- Optional for Web AOT: `dotnet workload install wasm-tools`

### Build the Solution
```bash
dotnet build TwistNTurn.sln
```

### Run Tests
```bash
dotnet test TwistNTurn.sln
```

### Run the Web Application (Blazor)
```bash
dotnet run --project TwistNTurn.Web
```
Then open `http://localhost:5000` (or the URL displayed in the terminal) in your browser.

### Publish for Production (with WebAssembly AOT)
```bash
dotnet publish TwistNTurn.Web/TwistNTurn.Web.csproj -c Release -o publish
```
The files in `publish/wwwroot` can be statically hosted anywhere (e.g. GitHub Pages, Cloudflare Pages, Netlify, Nginx, Apache).

## Game Controls & Features (Web App)

- **Left-Click / Drag edge**: Cycle state (Empty ➔ Line ➔ Cross ➔ Empty).
- **Right-Click / Drag edge**: Mark Excluded crosses (X).
- **Shift + Click**: Reverses cycling order for clicks.
- **Cell Shading**: When enabled, clicking a cell toggles background shading (Yellow / Gray) to help track interior vs. exterior loop parity.
- **Restrict Invalid Option**: Restricts manual edge marks from causing immediate local contradictions (vertex degree > 2).
- **Checkpoints**:
  - **Fix** (`Ctrl+F`): Lock and color current moves green as a known-good checkpoint. Fixed edges cannot be accidentally modified.
  - **Unfix** (`Ctrl+U`): Unlock fixed edges so they can be changed again.
  - **Revert** (`Ctrl+R`): Discard subsequent moves and roll back to the saved checkpoint.
- **History**: **Undo** (`Ctrl+Z`) and **Redo** (`Ctrl+Y`).
- **Touch & Gestures**: Two-finger pinch to zoom and drag to pan; single-finger touch draws loop lines.
- **Share**: Copy a compact URL link containing the complete puzzle state to share with friends.

