import type { GraphEdge, GraphNode } from "../types/cytocarto";

type PositionedNode = GraphNode & { x: number; y: number };

export const NETWORK_ARROW_LENGTH = 10;
export const NETWORK_ARROW_SIZE = 11;
const NETWORK_ARROW_TIP_GAP = 3;

export type EdgeRoute = {
  sourceOffset?: number;
  targetOffset?: number;
};

const NODE_RADIUS: Record<GraphNode["kind"], number> = {
  analyte: 19,
  gene: 22,
  program: 31,
  cell: 32,
};

export function centeredOffsets(count: number): number[] {
  if (count <= 1) return [0];
  const spacing = Math.min(8, 24 / (count - 1));
  return Array.from({ length: count }, (_, index) => (index - (count - 1) / 2) * spacing);
}

function boundaryPoint(node: PositionedNode, unitX: number, unitY: number, offset: number, forward: boolean, gap = 0) {
  const radius = NODE_RADIUS[node.kind] + gap;
  const limitedOffset = Math.max(-radius * 0.75, Math.min(radius * 0.75, offset));
  const radialDistance = Math.sqrt(Math.max(0, radius * radius - limitedOffset * limitedOffset));
  const direction = forward ? 1 : -1;
  return {
    x: node.x + direction * unitX * radialDistance - unitY * limitedOffset,
    y: node.y + direction * unitY * radialDistance + unitX * limitedOffset,
  };
}

function horizontalBoundaryPoint(node: PositionedNode, direction: number, offset: number, gap = 0) {
  const radius = NODE_RADIUS[node.kind] + gap;
  const limitedOffset = Math.max(-radius * 0.75, Math.min(radius * 0.75, offset));
  return {
    x: node.x + direction * Math.sqrt(Math.max(0, radius * radius - limitedOffset * limitedOffset)),
    y: node.y + limitedOffset,
  };
}

export function edgePath(edge: GraphEdge, source: PositionedNode, target: PositionedNode, route: EdgeRoute = {}): string {
  const distance = Math.hypot(target.x - source.x, target.y - source.y);
  const unitX = (target.x - source.x) / (distance || 1);
  const unitY = (target.y - source.y) / (distance || 1);
  const targetGap = edge.kind === "cell_context" ? 0 : NETWORK_ARROW_LENGTH + NETWORK_ARROW_TIP_GAP;
  const horizontalDirection = Math.sign(target.x - source.x);
  const start = horizontalDirection && edge.kind !== "program_program"
    ? horizontalBoundaryPoint(source, horizontalDirection, route.sourceOffset ?? 0)
    : boundaryPoint(source, unitX, unitY, route.sourceOffset ?? 0, true);
  const end = horizontalDirection && edge.kind !== "program_program"
    ? horizontalBoundaryPoint(target, -horizontalDirection, route.targetOffset ?? 0, targetGap)
    : boundaryPoint(target, unitX, unitY, route.targetOffset ?? 0, false, targetGap);
  if (edge.kind === "program_program") {
    const bend = source.y < target.y ? 70 : -70;
    return `M ${start.x} ${start.y} C ${start.x + bend} ${start.y}, ${end.x + bend} ${end.y}, ${end.x} ${end.y}`;
  }
  const middle = (start.x + end.x) / 2;
  return `M ${start.x} ${start.y} C ${middle} ${start.y}, ${middle} ${end.y}, ${end.x} ${end.y}`;
}
