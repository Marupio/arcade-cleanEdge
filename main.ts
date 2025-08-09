/*** MIT LICENSE
Copyright (c) 2022 torcado
Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:
The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
***/

//% color=#C040C0 icon="\uf074" block="CleanEdge Rotate"
namespace cleanedgerotate {
    // ------- Options -------
    export interface Options {
        // how similar two colours must be to be considered same (Euclidean in RGBA [0..1])
        similarThreshold?: number;   // default 0
        // perceived stroke width (see shader notes)
        lineWidth?: number;          // default 1.0
        // slope/cleanup toggles (speed vs quality)
        enableSlope?: boolean;       // default true
        enableCleanup?: boolean;     // default false (rotation-only use case)
        // colour priority reference (palette index). 0 is transparent.
        highestColorIndex?: number;  // default 1 (white in default palette)
        // centre override (default: image centre)
        cx?: number;
        cy?: number;
        // output size (defaults to source size)
        outW?: number;
        outH?: number;
    }

    // -------- Palette: Arcade’s default 16-colour palette (index -> RGB) --------
    // 0 == transparent; treat as alpha=0.
    const PALETTE: number[] = [
        0x000000, 0xFFFFFF, 0xFF2121, 0xFF93C4, 0xFF8135, 0xFFF609, 0x249CA3, 0x78DC52,
        0x003FAD, 0x87F2FF, 0x8E2EC4, 0xA4839F, 0x5C406C, 0xE5CDC4, 0x91463D, 0x000000
    ];

    // Convert palette index -> {r,g,b,a} in [0..1]
    function col4(idx: number): number[] {
        if (idx <= 0) return [0, 0, 0, 0]; // transparent
        const rgb = PALETTE[idx & 15];
        const r = ((rgb >> 16) & 0xFF) / 255;
        const g = ((rgb >> 8) & 0xFF) / 255;
        const b = (rgb & 0xFF) / 255;
        return [r, g, b, 1];
    }

    function distRGBA(a: number[], b: number[]): number {
        const dr = a[0]-b[0], dg = a[1]-b[1], db = a[2]-b[2], da = a[3]-b[3];
        return Math.sqrt(dr*dr + dg*dg + db*db + da*da);
    }

    function similar(a: number[], b: number[], thr: number): boolean {
        // transparent counts as same as transparent
        if (a[3] === 0 && b[3] === 0) return true;
        return distRGBA(a, b) <= thr;
    }

    function higher(thisCol: number[], other: number[], highestRGB: number[]): boolean {
        // if similar, not higher
        // alpha priority; else distance to highestColor (RGB) smaller wins
        if (thisCol[3] === 0 && other[3] === 0) return false;
        if (thisCol[3] === other[3]) {
            const d1 = Math.sqrt(
                (thisCol[0]-highestRGB[0])**2 +
                (thisCol[1]-highestRGB[1])**2 +
                (thisCol[2]-highestRGB[2])**2
            );
            const d2 = Math.sqrt(
                (other[0]-highestRGB[0])**2 +
                (other[1]-highestRGB[1])**2 +
                (other[2]-highestRGB[2])**2
            );
            return d1 < d2;
        }
        return thisCol[3] > other[3];
    }

    // signed distance from point to line through pt1->pt2; sign by side of `dir`
    function distToLine(px: number, py: number, x1: number, y1: number, x2: number, y2: number, dirx: number, diry: number): number {
        const lx = x2 - x1, ly = y2 - y1;
        const perpx = ly, perpy = -lx;
        const dx = x1 - px, dy = y1 - py;
        const sign = (perpx*dirx + perpy*diry) > 0 ? 1 : -1;
        const plen = Math.sqrt(perpx*perpx + perpy*perpy) || 1;
        // normalize perp, then dot with dir to point
        return sign * ((perpx/plen)*dx + (perpy/plen)*dy);
    }

    function clamp(v: number, lo: number, hi: number) { return v < lo ? lo : (v > hi ? hi : v); }

    // Sample palette colour as RGBA vec at integer coordinate (out-of-bounds -> transparent)
    function sample(img: Image, x: number, y: number): number[] {
        if (x < 0 || y < 0 || x >= img.width || y >= img.height) return [0,0,0,0];
        return col4(img.getPixel(x, y));
    }

    // Helper: fetch 5x5 neighbourhood around (ix,iy), arranged to match shader variable names
    // ub u  uf uff
    //  b c  f  ff
    // db d  df dff
    //ddb dd ddf
    function neigh(img: Image, ix: number, iy: number) {
        const g = (dx: number, dy: number) => sample(img, ix+dx, iy+dy);
        return {
            ub: g(0,-2), u: g(0,-1), uf: g(1,-1), uff: g(2,-1),
            b: g(-1,0), c: g(0,0), f: g(1,0), ff: g(2,0),
            db: g(-1,1), d: g(0,1), df: g(1,1), dff: g(2,1),
            ddb: g(-1,2), dd: g(0,2), ddf: g(1,2)
        };
    }

    function defnum(v: number, d: number): number {
        return (v === undefined || v === null) ? d : v;
    }
    function defbool(v: boolean, d: boolean): boolean {
        return (v === undefined || v === null) ? d : v;
    }
    function sqr(x: number): number { return x * x; } // safer than x**2 in Arcade


    // Core of the shader: given the local neighbourhood and a sub-pixel point within the current texel,
    // decide whether to take c, f, d, etc., using the edge rules (45° + optional 2:1 slopes).
    function sliceDecide(
        px: number, py: number, // subpixel point in [0..1]^2 within the centre cell
        dirMainX: number, dirMainY: number,
        dirPtX: number, dirPtY: number,
        nb: any,
        options: Options,
        highestRGB: number[]
    ): number[] {
        const thr = defnum(options.similarThreshold, 0.0);
        const enableSlope = defbool(options.enableSlope, true);
        const enableCleanup = defbool(options.enableCleanup, false);
        const LW = defnum(options.lineWidth, 1.0);

        const minW = enableSlope ? 0.45 : 0.0;
        const maxW = enableSlope ? 1.142 : 1.4;
        const lineW = clamp(LW, minW, maxW);

        // flip point like shader
        const fx = dirMainX * (px - 0.5) + 0.5;
        const fy = dirMainY * (py - 0.5) + 0.5;

        // Small helpers
        const cd = (A: number[], B: number[]) => distRGBA(A, B);
        const sim = (A: number[], B: number[]) => (A[3]===0 && B[3]===0) || cd(A,B) <= thr;

        // Common vars
        const {ub,u,uf,uff,b,c,f,ff,db,d,df,dff,ddb,dd,ddf} = nb;
        const centerX = 0.5, centerY = 0.5;

        // Edge detection
        const distAgainst = 4.0*cd(f,d) + cd(uf,c) + cd(c,db) + cd(ff,df) + cd(df,dd);
        const distTowards = 4.0*cd(c,df) + cd(u,f) + cd(f,dff) + cd(b,d) + cd(d,ddf);
        let shouldSlice = (distAgainst < distTowards) || ((distAgainst < distTowards + 0.001) && !higher(c,f,highestRGB));
        if (sim(f,d) && sim(b,u) && sim(uf,df) && sim(db,ub) && !sim(c,f)) {
            // checkerboard edge case
            shouldSlice = false;
        }
        if (!shouldSlice) return [-1,-1,-1,-1];

        let dist = 1.0;
        let flip = false;

        // Helper for distToLine with direction
        const dline = (ax:number, ay:number, bx:number, by:number, sgn:number) =>
            distToLine(fx, fy, centerX+ax*dirPtX, centerY+ay*dirPtY, centerX+bx*dirPtX, centerY+by*dirPtY, sgn*dirPtX, sgn*dirPtY);

        // ---- Cases (mirroring shader structure) ----
        // To keep this readable and performant on Arcade, we implement the main three branches:
        // (1) 2:1 shallow/steep (if enableSlope), (2) 45° diagonal, (3) far-corner slope helpers (if enableSlope).
        // The numeric constants follow the shader.

        // Utility to finalize a choice between two candidates by distance to c (tie-break like shader)
        function chooseByProximity(A: number[], B: number[]): number[] {
            return cd(c, A) <= cd(c, B) ? A : B;
        }

        if (enableSlope && sim(f,d) && !sim(f,b) && !sim(uf,db)) {
            // 45° path will handle; fall through.
        }

        // ---------- 45° diagonal ----------
        if (sim(f,d)) {
            if (sim(c,df) && higher(c,f,highestRGB)) {
                if (!sim(c,dd) && !sim(c,ff)) flip = true;
            } else {
                if (higher(c,f,highestRGB)) flip = true;
                if (!sim(c,b) && sim(b,f) && sim(b,d) && sim(b,u)) flip = true;
            }

            // single-pixel 2:1 slope special-case
            if (((sim(f,db) && sim(u,f) && sim(f,df)) || (sim(uf,d) && sim(b,d) && sim(d,df))) && !sim(c,df))
                flip = true;

            if (flip) {
                dist = lineW - dline( 1.0,-1.0, -1.0, 1.0, -1);
            } else {
                dist = dline(1.0,0.0, 0.0,1.0, +1);
            }

            // optional cleanup (reduced to minimal influence for perf)
            dist -= (lineW/2.0);
            return dist <= 0.0 ? chooseByProximity(f, d) : [-1,-1,-1,-1];
        }

        // ---------- Optional 2:1 slopes (shallow & steep) ----------
        if (enableSlope) {
            // Lower shallow 2:1: similar(f,d,db) && !similar(f,d,b) && !similar(uf,db)
            if (sim(f,d) && sim(d,db) && !(sim(f,d) && sim(d,b)) && !sim(uf,db)) {
                if (!(sim(c,df) && higher(c,f,highestRGB))) {
                    if (higher(c,f,highestRGB)) flip = true;
                    if (sim(u,f) && !sim(c,df) && !higher(c,u,highestRGB)) flip = true;
                }
                if (flip) dist = lineW - dline(1.5,-1.0, -0.5,0.0, -1);
                else      dist = dline(1.5, 0.0, -0.5,1.0, +1);

                if (enableCleanup && !flip && sim(c,uf)) {
                    const d2 = dline(2.0,-1.0, 0.0,1.0, +1);
                    dist = Math.min(dist, d2);
                }
                dist -= (lineW/2.0);
                return dist <= 0.0 ? chooseByProximity(f, d) : [-1,-1,-1,-1];
            }

            // Forward steep 2:1: similar(uf,f,d) && !similar(u,f,d) && !similar(uf,db)
            if (sim(uf,f) && sim(f,d) && !(sim(u,f) && sim(f,d)) && !sim(uf,db)) {
                if (!(sim(c,df) && higher(c,d,highestRGB))) {
                    if (higher(c,d,highestRGB)) flip = true;
                    if (sim(b,d) && !sim(c,df) && !higher(c,d,highestRGB)) flip = true;
                }
                if (flip) dist = lineW - dline(0.0,-0.5, -1.0,1.5, -1);
                else      dist = dline(1.0,-0.5,  0.0,1.5, +1);

                if (enableCleanup && !flip && sim(c,db)) {
                    const d2 = dline(1.0,0.0, -1.0,2.0, +1);
                    dist = Math.min(dist, d2);
                }
                dist -= (lineW/2.0);
                return dist <= 0.0 ? chooseByProximity(f, d) : [-1,-1,-1,-1];
            }
        }

        // Fallback: no slice → no change
        return [-1,-1,-1,-1];
    }

    function indexOfColor(col: number[]): number {
        // Pick the nearest palette index for the chosen RGBA (ignoring alpha==0 -> 0)
        if (col[3] === 0) return 0;
        let best = 1, bestD = 1e9;
        for (let i = 1; i < 16; i++) {
            const c = col4(i);
            const d = (c[0]-col[0])**2 + (c[1]-col[1])**2 + (c[2]-col[2])**2;
            if (d < bestD) { bestD = d; best = i; }
        }
        return best;
    }

    // -------- Public API --------

    /**
     * Rotate src by angle (degrees), preserving crisp edges (cleanEdge-inspired).
     * @param src source image
     * @param angleDeg rotation in degrees (clockwise positive)
     * @param opts optional tuning parameters
     */
    //% blockId=cleanedge_rotate block="rotate (clean edge) %src by %angleDeg (°)"
    //% src.shadow=variables_get angleDeg.defl=30 weight=100
    export function rotateCleanEdge(src: Image, angleDeg: number, opts?: Options): Image {
        if (!src) return null;
        const options = opts || {};
        const outW = defnum(options.outW, src.width);
        const outH = defnum(options.outH, src.height);
        const cx = defnum(options.cx, (src.width - 1) / 2);
        const cy = defnum(options.cy, (src.height - 1) / 2);

        // highest colour reference
        const hiIdx = defnum(options.highestColorIndex, 1);
        const hiRGB = col4(hiIdx); // use RGB from palette (alpha=1)

        const dst = image.create(outW, outH);
        dst.fill(0);

        const a = angleDeg * Math.PI / 180;
        const cosA = Math.cos(a), sinA = Math.sin(a);

        for (let y = 0; y < outH; y++) {
            const dy = y - cy;
            for (let x = 0; x < outW; x++) {
                const dx = x - cx;

                // inverse rotate
                const sx =  cosA * dx + sinA * dy + cx;
                const sy = -sinA * dx + cosA * dy + cy;

                // integer centre cell
                const ix = Math.floor(sx);
                const iy = Math.floor(sy);
                const fx = sx - ix; // subpixel in [0..1)
                const fy = sy - iy;

                // pull neighbourhood
                const nb = neigh(src, ix, iy);

                // Try four principal directions (like shader’s quadrant logic)
                // Using unit vectors for main/point direction permutations:
                const candidates: number[][] = [];

                // down-forward basis
                const c1 = sliceDecide(fx, fy, 1, 1, 1, 1, nb, options, hiRGB);
                if (c1[0] >= 0) candidates.push(c1);

                // down-back
                const c2 = sliceDecide(1-fx, fy, -1, 1, -1, 1, nb, options, hiRGB);
                if (c2[0] >= 0) candidates.push(c2);

                // up-forward
                const c3 = sliceDecide(fx, 1-fy, 1, -1, 1, -1, nb, options, hiRGB);
                if (c3[0] >= 0) candidates.push(c3);

                // up-back
                const c4 = sliceDecide(1-fx, 1-fy, -1, -1, -1, -1, nb, options, hiRGB);
                if (c4[0] >= 0) candidates.push(c4);

                let outCol: number[];
                if (candidates.length) {
                    // start with first element, then compare to rest
                    let best = candidates[0];
                    for (let i = 1; i < candidates.length; i++) {
                        const cur = candidates[i];
                        if (higher(cur, best, hiRGB)) best = cur;
                    }
                    outCol = best;
                } else {
                    // No slice override → nearest neighbour
                    const sxr = Math.round(sx);
                    const syr = Math.round(sy);
                    outCol = sample(src, sxr, syr);
                }

                dst.setPixel(x, y, indexOfColor(outCol));
            }
        }
        return dst;
    }
}
