# arcade-cleanedge-rotate

Edge-preserving rotation for pixel art in MakeCode Arcade, inspired by torcado's cleanEdge shader (MIT).
Exposes a `rotateCleanEdge()` API and blocks.

### Options
- `similarThreshold` (0..~0.25): how close colours must be to count as “same”
- `lineWidth` (~0.45–1.14 with slope, 0–1.4 otherwise): perceived stroke width
- `enableSlope`: handle 2:1 slants (faster off)
- `enableCleanup`: clean small slope transitions (faster off)
- `highestColorIndex`: priority reference (default: white)

### Example
```ts
const rotated = cleanedgerotate.rotateCleanEdge(myImg, 30)
sprites.create(rotated)
