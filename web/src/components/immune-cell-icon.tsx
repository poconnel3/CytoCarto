type CellIconKind =
  | "activated-macrophage"
  | "b-cell"
  | "basophil"
  | "cancer"
  | "dendritic"
  | "eosinophil"
  | "epithelial"
  | "generic"
  | "macrophage"
  | "mast-cell"
  | "monocyte"
  | "natural-killer"
  | "neutrophil"
  | "stem-cell"
  | "t-cell";

type Crop = {
  src: string;
  imageWidth: number;
  imageHeight: number;
  x: number;
  y: number;
  size: number;
};

const PRIMARY = "/immune-cell-symbols-primary.png";
const SECONDARY = "/immune-cell-symbols-secondary.png";

const CROPS: Record<CellIconKind, Crop> = {
  "activated-macrophage": { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 195, y: 6, size: 148 },
  "b-cell": { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 359, y: 246, size: 148 },
  basophil: { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 358, y: 5, size: 148 },
  cancer: { src: SECONDARY, imageWidth: 1108, imageHeight: 448, x: 21, y: 14, size: 142 },
  dendritic: { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 519, y: 0, size: 152 },
  eosinophil: { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 682, y: 4, size: 150 },
  epithelial: { src: SECONDARY, imageWidth: 1108, imageHeight: 448, x: 974, y: 10, size: 132 },
  generic: { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 1010, y: 246, size: 150 },
  macrophage: { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 27, y: 9, size: 148 },
  "mast-cell": { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 842, y: 242, size: 156 },
  monocyte: { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 844, y: 6, size: 148 },
  "natural-killer": { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 1007, y: 4, size: 148 },
  neutrophil: { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 682, y: 245, size: 150 },
  "stem-cell": { src: SECONDARY, imageWidth: 1108, imageHeight: 448, x: 809, y: 12, size: 144 },
  "t-cell": { src: PRIMARY, imageWidth: 1168, imageHeight: 512, x: 28, y: 247, size: 148 },
};

export function immuneCellKind(cellType: string): CellIconKind {
  const value = cellType.toLowerCase();
  if (value.includes("activated macrophage") || value.includes("inflammatory macrophage")) return "activated-macrophage";
  if (value.includes("macrophage")) return "macrophage";
  if (value.includes("monocyte") || value.includes("mono")) return "monocyte";
  if (value.includes("dendritic") || value.includes("cdc") || value.includes("pdc") || value === "dc") return "dendritic";
  if (value.includes("natural killer") || /(^|\W)nk($|\W)/.test(value)) return "natural-killer";
  if (value.includes("neutrophil")) return "neutrophil";
  if (value.includes("eosinophil")) return "eosinophil";
  if (value.includes("basophil")) return "basophil";
  if (value.includes("mast")) return "mast-cell";
  if (value.includes("cancer") || value.includes("tumor")) return "cancer";
  if (value.includes("epithelial")) return "epithelial";
  if (value.includes("stem") || value.includes("hsc")) return "stem-cell";
  if (value.includes("b cell") || value.includes("b-cell") || value.includes("plasma")) return "b-cell";
  if (value.includes("t cell") || value.includes("t-cell") || value.includes("cd4") || value.includes("cd8") || value.includes("treg")) return "t-cell";
  return "generic";
}

export function ImmuneCellIcon({
  cellType,
  size = 42,
  x,
  y,
  className = "",
}: {
  cellType: string;
  size?: number;
  x?: number;
  y?: number;
  className?: string;
}) {
  const crop = CROPS[immuneCellKind(cellType)];
  return <svg
    className={`immune-cell-symbol ${className}`}
    x={x}
    y={y}
    width={size}
    height={size}
    viewBox={`${crop.x} ${crop.y} ${crop.size} ${crop.size}`}
    aria-hidden="true"
    focusable="false"
  >
    <image href={crop.src} width={crop.imageWidth} height={crop.imageHeight} />
  </svg>;
}
