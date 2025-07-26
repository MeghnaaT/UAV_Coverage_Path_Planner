# Add this at the top with your imports
from shapely.geometry import Point
from itertools import permutations
import matplotlib.pyplot as plt
import re  
from xml.dom import minidom
from shapely.geometry import Polygon
from shapely.ops import triangulate, unary_union
from shapely.geometry import LineString, MultiLineString
import simplekml
import matplotlib.pyplot as plt
import numpy as np


# Set the values as per Hawkeye Thumb 4K
from math import radians, tan, sin

def get_path_spacing(altitude, horizontal_fov=168.5, vertical_fov=159.76, camera_pitch=90, overlap=0.2):
    horizontal_fov = radians(horizontal_fov)
    vertical_fov = radians(vertical_fov)
    camera_pitch = radians(camera_pitch)

    spacing = (2 * altitude * tan(horizontal_fov / 2)) / (sin(camera_pitch + vertical_fov / 2))
    spacing = spacing - (spacing * overlap)
    return spacing

def extract_polygon_from_kml(kml_file_path):
    """Extract polygon from KML file"""
    kml_doc = minidom.parse(kml_file_path)
    coord_elements = kml_doc.getElementsByTagName("coordinates")

    for coords in coord_elements:
        coord_text = coords.firstChild.nodeValue.strip()
        coord_pairs = re.findall(r"([-.\d]+),([-.\d]+)", coord_text)
        latlon = [(float(lon), float(lat)) for lon, lat in coord_pairs]
        if len(latlon) > 2:
            polygon = Polygon(latlon)
            if polygon.is_valid:
                return polygon
    raise ValueError("No valid polygon found in KML.")


def is_convex(polygon, tolerance=1e-10):
    """Check if polygon is convex using cross product method"""
    if not polygon.is_valid:
        return False

    coords = list(polygon.exterior.coords)[:-1]  # Remove duplicate closing point
    if len(coords) < 3:
        return True

    def cross_product_z(a, b, c):
        """Calculate z-component of cross product for 2D points"""
        ab = (b[0] - a[0], b[1] - a[1])
        bc = (c[0] - b[0], c[1] - b[1])
        return ab[0] * bc[1] - ab[1] * bc[0]

    # Check all consecutive triplets of vertices
    signs = []
    for i in range(len(coords)):
        a = coords[i]
        b = coords[(i + 1) % len(coords)]
        c = coords[(i + 2) % len(coords)]
        cp = cross_product_z(a, b, c)
        if abs(cp) > tolerance:  # Ignore near-zero cross products
            signs.append(cp > 0)

    # All non-zero cross products should have same sign for convex polygon
    return len(set(signs)) <= 1


def hertel_mehlhorn_decomposition(polygon):
    """
    Hertel-Mehlhorn algorithm for convex decomposition
    This ensures complete coverage by working with the polygon structure directly
    """
    if is_convex(polygon):
        return [polygon]

    # Start with triangulation but use a different approach
    triangles = list(triangulate(polygon))

    # Find triangles that are completely or mostly inside the polygon
    valid_triangles = []
    for tri in triangles:
        intersection = polygon.intersection(tri)
        if intersection.area > 0.01 * tri.area:  # At least 1% overlap
            if intersection.geom_type == "Polygon":
                valid_triangles.append(intersection)
            else:
                for geom in intersection.geoms:
                    if geom.geom_type == "Polygon" and geom.area > 1e-10:
                        valid_triangles.append(geom)

    # Remove duplicates and very small triangles
    filtered_triangles = []
    min_area = polygon.area * 1e-6

    for tri in valid_triangles:
        if tri.area > min_area:
            # Check if this triangle is already covered by existing triangles
            is_duplicate = False
            for existing in filtered_triangles:
                if (
                    abs(tri.area - existing.area) < min_area
                    and tri.centroid.distance(existing.centroid) < 1e-8
                ):
                    is_duplicate = True
                    break
            if not is_duplicate:
                filtered_triangles.append(tri)

    # Add any missing areas by creating additional triangles
    union_triangles = (
        unary_union(filtered_triangles) if filtered_triangles else Polygon()
    )
    missing_area = polygon.difference(union_triangles)

    if missing_area.area > min_area:
        if missing_area.geom_type == "Polygon":
            # Triangulate the missing area
            try:
                missing_triangles = list(triangulate(missing_area))
                for tri in missing_triangles:
                    if missing_area.contains(tri.centroid):
                        filtered_triangles.append(tri)
            except:
                # If triangulation fails, add the missing area as is
                filtered_triangles.append(missing_area)
        elif hasattr(missing_area, "geoms"):
            for geom in missing_area.geoms:
                if geom.geom_type == "Polygon" and geom.area > min_area:
                    try:
                        missing_triangles = list(triangulate(geom))
                        for tri in missing_triangles:
                            if geom.contains(tri.centroid):
                                filtered_triangles.append(tri)
                    except:
                        filtered_triangles.append(geom)

    return merge_triangles_to_convex(filtered_triangles, polygon)


def merge_triangles_to_convex(triangles, original_polygon):
    """Merge triangles into larger convex regions"""
    if not triangles:
        return [original_polygon] if is_convex(original_polygon) else []

    regions = triangles.copy()

    # Multiple passes of merging
    for iter in range(10):
        new_regions = []
        used = set()

        # Sort by area (largest first)
        regions.sort(key=lambda x: x.area, reverse=True)

        for i, region1 in enumerate(regions):
            if i in used:
                continue

            best_merge = region1
            best_partner = None

            # Try to merge with other regions
            for j in range(i + 1, len(regions)):
                if j in used:
                    continue

                region2 = regions[j]

                # Check if regions are adjacent (share a boundary)
                if region1.touches(region2) or region1.intersects(region2):
                    try:
                        merged = unary_union([region1, region2])
                        if (
                            merged.geom_type == "Polygon"
                            and merged.is_valid
                            and is_convex(merged)
                            and merged.area > best_merge.area
                        ):
                            best_merge = merged
                            best_partner = j
                    except:
                        continue

            if best_partner is not None:
                used.add(i)
                used.add(best_partner)
                new_regions.append(best_merge)
            else:
                used.add(i)
                new_regions.append(region1)

        if len(new_regions) >= len(regions):
            break  # No improvement

        regions = new_regions

    return regions


def optimize_convex_decomposition(polygon, buffer_space_degrees, spacing_degrees):

    """Main function for convex decomposition with complete coverage

    Args:
        polygon: The polygon to decompose (possibly buffered)
        original_polygon: The original polygon before buffering (optional, for visualization)
    """
    original_polygon = polygon
    
    polygon = polygon.buffer(
    -buffer_space_degrees,
    quad_segs=2,
    cap_style="round",
    join_style="mitre",
    mitre_limit=1,
)
  
    if is_convex(polygon):
        print("Polygon is already convex!")
        # Visualization for convex polygon
        fig, ax = plt.subplots(1, 1, figsize=(16, 8))
        ax.set_title("Convex Polygon (No Decomposition Needed)")
        if original_polygon is not None:
            x0, y0 = original_polygon.exterior.xy
            ax.plot(x0, y0, "b--", linewidth=2, label="Original Polygon (pre-buffer)")
        x, y = polygon.exterior.xy
        ax.plot(x, y, "k-", linewidth=2, label="Buffered Polygon")
        ax.fill(x, y, color="lightgreen", alpha=0.6, edgecolor="green", linewidth=2)
        ax.set_aspect("equal")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        return [polygon]

    # Use Hertel-Mehlhorn based approach
    convex_regions = hertel_mehlhorn_decomposition(polygon)

    # Verify and fix coverage
    total_union = unary_union(convex_regions) if convex_regions else Polygon()
    coverage_ratio = total_union.area / polygon.area if polygon.area > 0 else 0

    print(f"Initial coverage: {coverage_ratio * 100:.2f}%")

    # If coverage is incomplete, add missing areas
    if coverage_ratio < 0.99:
        missing = polygon.difference(total_union)
        if missing.area > polygon.area * 1e-6:
            if missing.geom_type == "Polygon":
                if is_convex(missing):
                    convex_regions.append(missing)
                else:
                    # Recursively decompose missing area
                    missing_regions = hertel_mehlhorn_decomposition(missing)
                    convex_regions.extend(missing_regions)
            elif hasattr(missing, "geoms"):
                for geom in missing.geoms:
                    if geom.geom_type == "Polygon" and geom.area > polygon.area * 1e-6:
                        if is_convex(geom):
                            convex_regions.append(geom)
                        else:
                            missing_regions = hertel_mehlhorn_decomposition(geom)
                            convex_regions.extend(missing_regions)

    # Final verification and cleanup
    final_regions = []
    min_area = polygon.area * 1e-6

    for region in convex_regions:
        if region.area > min_area and region.is_valid:
            final_regions.append(region)

    # --- Merge small polygons with neighbors ---
    merge_min_area = 0.0000005
    merged = True
    while merged:
        merged = False
        small_regions = [r for r in final_regions if r.area < merge_min_area]
        if not small_regions:
            break
        for small in small_regions:
            # Find neighbors (touching polygons)
            neighbors = [r for r in final_regions if r != small and r.touches(small)]
            if neighbors:
                # Pick neighbor with largest shared boundary
                neighbor = max(neighbors, key=lambda n: small.intersection(n).length)
                new_region = unary_union([small, neighbor])
                final_regions.remove(small)
                final_regions.remove(neighbor)
                final_regions.append(new_region)
                merged = True
                break  # Restart loop since list changed

    # Visualization
    fig, ax = plt.subplots(1, 1, figsize=(16, 8))

    ax.set_title(f"Convex Decomposition ({len(final_regions)} regions)")
    if original_polygon is not None:
        x0, y0 = original_polygon.exterior.xy
        ax.plot(x0, y0, "b--", linewidth=2, label="Original Polygon (pre-buffer)")
    x, y = polygon.exterior.xy
    ax.plot(x, y, "k-", linewidth=2, label="Buffered Polygon")

    colors = plt.cm.Set3(np.linspace(0, 1, max(len(final_regions), 1)))
    for i, region in enumerate(final_regions):
        x, y = region.exterior.xy
        color = colors[i % len(colors)]
        ax.fill(x, y, color=color, alpha=0.6, edgecolor="red", linewidth=1)

    # --- Draw lawnmower path lines ---
    paths = generate_lawnmower_path_for_cells(final_regions, spacing_degrees)
    for path in paths:
        if hasattr(path, 'coords'):
            x, y = zip(*path.coords)
            ax.plot(x, y, color='black', linewidth=1.2, alpha=0.8, linestyle='--', label=None)

    ax.set_aspect("equal")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    # Final verification
    final_union = unary_union(final_regions) if final_regions else Polygon()
    total_area = final_union.area
    coverage = total_area / polygon.area * 100 if polygon.area > 0 else 0

    print(f"\nFinal Decomposition Results:")
    print(f"Number of convex regions: {len(final_regions)}")
    print(f"Total area coverage: {coverage:.2f}%")
    print(f"Original area: {polygon.area:.6f}")
    print(f"Decomposed area: {total_area:.6f}")

    for i, region in enumerate(final_regions):
        print(f"Region {i+1}: Area = {region.area:.6f}, Convex = {is_convex(region)}")

    return final_regions


def export_convex_cells_to_kml(convex_cells, original_polygon, output_path):
    """Export convex decomposition to KML file"""
    kml = simplekml.Kml()

    # Original polygon in blue
    pol = kml.newpolygon(name="Original Polygon")
    pol.outerboundaryis = [(lon, lat) for lon, lat in original_polygon.exterior.coords]
    pol.style.polystyle.color = simplekml.Color.changealphaint(
        100, simplekml.Color.blue
    )
    pol.style.linestyle.color = simplekml.Color.blue
    pol.style.linestyle.width = 2

    # Convex regions in different colors
    colors = [
        simplekml.Color.red,
        simplekml.Color.green,
        simplekml.Color.yellow,
        simplekml.Color.purple,
        simplekml.Color.orange,
        simplekml.Color.pink,
    ]

    for i, cell in enumerate(convex_cells):
        coords = [(lon, lat) for lon, lat in cell.exterior.coords]
        p = kml.newpolygon(name=f"Convex Region {i+1}")
        p.outerboundaryis = coords
        color = colors[i % len(colors)]
        p.style.polystyle.color = simplekml.Color.changealphaint(120, color)
        p.style.linestyle.color = color
        p.style.linestyle.width = 2

    kml.save(output_path)
    print(f"Convex decomposition saved to: {output_path}")

def generate_lawnmower_path_for_cells(convex_cells, spacing_deg):
    """Generate lawnmower path in each convex cell, starting with the first line passing through the longest edge."""
    all_paths = []

    for idx, cell in enumerate(convex_cells):
        coords = list(cell.exterior.coords)
        # 1. Find the longest edge
        max_len = 0
        longest_edge = (coords[0], coords[1])
        for i in range(len(coords) - 1):
            p1 = coords[i]
            p2 = coords[i + 1]
            length = ((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2) ** 0.5
            if length > max_len:
                max_len = length
                longest_edge = (p1, p2)

        (x1, y1), (x2, y2) = longest_edge
        direction = ((x2 - x1) / max_len, (y2 - y1) / max_len)
        perp_dir = (-direction[1], direction[0])

        # Project all points onto the perpendicular direction to get bounds
        projections = [p[0]*perp_dir[0] + p[1]*perp_dir[1] for p in coords]
        min_proj = min(projections)
        max_proj = max(projections)

        sweep_lines = []

        # --- First sweep line: exactly colinear with the longest edge ---
        # Extend the longest edge far in both directions
        dx, dy = direction
        # Use the edge's midpoint as base, but extend from both endpoints
        line = LineString([
            (x1 - 10000*dx, y1 - 10000*dy),
            (x2 + 10000*dx, y2 + 10000*dy)
        ])
        clipped = cell.buffer(0.00001).intersection(line)
        if not clipped.is_empty:
            if isinstance(clipped, LineString):
                sweep_lines.append(clipped)
            elif isinstance(clipped, MultiLineString):
                sweep_lines.extend([seg for seg in clipped.geoms if isinstance(seg, LineString)])

        # --- Now generate all other sweep lines, offset from the first ---
        # The anchor projection is the projection of the first endpoint
        anchor_proj = x1 * perp_dir[0] + y1 * perp_dir[1]

        # Positive direction
        pos = anchor_proj + spacing_deg
        while pos <= max_proj + 1e-10:
            base_x = x1 + perp_dir[0] * (pos - anchor_proj)
            base_y = y1 + perp_dir[1] * (pos - anchor_proj)
            line = LineString([
                (base_x - 10000*dx, base_y - 10000*dy),
                (base_x + 10000*dx, base_y + 10000*dy)
            ])
            clipped = cell.intersection(line)
            if not clipped.is_empty:
                if isinstance(clipped, LineString):
                    sweep_lines.append(clipped)
                elif isinstance(clipped, MultiLineString):
                    sweep_lines.extend([seg for seg in clipped.geoms if isinstance(seg, LineString)])
            pos += spacing_deg

        # Negative direction
        pos = anchor_proj - spacing_deg
        while pos >= min_proj - 1e-10:
            base_x = x1 + perp_dir[0] * (pos - anchor_proj)
            base_y = y1 + perp_dir[1] * (pos - anchor_proj)
            line = LineString([
                (base_x - 10000*dx, base_y - 10000*dy),
                (base_x + 10000*dx, base_y + 10000*dy)
            ])
            clipped = cell.intersection(line)
            if not clipped.is_empty:
                if isinstance(clipped, LineString):
                    sweep_lines.append(clipped)
                elif isinstance(clipped, MultiLineString):
                    sweep_lines.extend([seg for seg in clipped.geoms if isinstance(seg, LineString)])
            pos -= spacing_deg

        # Sort sweep lines by their projection to maintain order
        def line_proj(line):
            mx = (line.coords[0][0] + line.coords[-1][0]) / 2
            my = (line.coords[0][1] + line.coords[-1][1]) / 2
            return mx * perp_dir[0] + my * perp_dir[1]
        sweep_lines = sorted(sweep_lines, key=line_proj)

        # Zigzag order
        ordered_coords = []
        flip = False
        for segment in sweep_lines:
            seg_coords = list(segment.coords)
            if flip:
                seg_coords = seg_coords[::-1]
            ordered_coords.extend(seg_coords)
            flip = not flip
        if ordered_coords:
            zigzag_line = LineString(ordered_coords)
            all_paths.append(zigzag_line)

    return all_paths

def export_lawnmower_path_to_kml(paths, convex_cells, output_file):
    kml = simplekml.Kml()

    # Draw convex cell boundaries
    for i, cell in enumerate(convex_cells):
        boundary = kml.newpolygon(name=f"Convex Cell {i+1}")
        coords = [(lon, lat) for lon, lat in cell.exterior.coords]
        boundary.outerboundaryis = coords
        boundary.style.polystyle.color = simplekml.Color.changealphaint(60, simplekml.Color.red)
        boundary.style.linestyle.color = simplekml.Color.red
        boundary.style.linestyle.width = 2

    # Draw lawnmower paths
    for i, path in enumerate(paths):
        coords = [(lon, lat) for lon, lat in path.coords]
        linestring = kml.newlinestring(name=f"Lawnmower Path {i+1}")
        linestring.coords = coords
        linestring.style.linestyle.color = simplekml.Color.green
        linestring.style.linestyle.width = 2

    kml.save(output_file)
    print(f"Lawnmower path with convex cells saved to: {output_file}")

# Add this at the top with your imports
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import nearest_points

# --- Helper: Distance ---
def euclidean(p1, p2):
    return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2) ** 0.5

# --- Helper: Choose direction ---
def choose_coverage_direction(current_pos, coverage_path):
    p0 = coverage_path[0]
    pn = coverage_path[-1]
    if euclidean(current_pos, p0) <= euclidean(current_pos, pn):
        return p0, pn, coverage_path  # forward
    else:
        return pn, p0, list(reversed(coverage_path))  # reverse

# --- Helper: A* transition planner constrained to polygon ---
def a_star_path(start, goal, polygon, step=0.5/111111):
    from heapq import heappush, heappop

    def neighbors(p):
        x, y = p
        return [
            (x + dx, y + dy)
            for dx in [-step, 0, step]
            for dy in [-step, 0, step]
            if (dx != 0 or dy != 0)
        ]

    def is_valid(p):
        return polygon.contains(Point(p))

    open_set = []
    came_from = {}
    g_score = {start: 0}
    f_score = {start: euclidean(start, goal)}
    heappush(open_set, (f_score[start], start))

    while open_set:
        _, current = heappop(open_set)
        if euclidean(current, goal) < step:
            # Reconstruct path
            path = [current]
            while current in came_from:
                current = came_from[current]
                path.append(current)
            return path[::-1]

        for neighbor in neighbors(current):
            if not is_valid(neighbor):
                continue
            tentative_g = g_score[current] + euclidean(current, neighbor)
            if neighbor not in g_score or tentative_g < g_score[neighbor]:
                came_from[neighbor] = current
                g_score[neighbor] = tentative_g
                f_score[neighbor] = tentative_g + euclidean(neighbor, goal)
                heappush(open_set, (f_score[neighbor], neighbor))

    return [start, goal]  # fallback if no path found

# --- Extract lawnmower path coordinates from LineString ---
def extract_path_coords(paths):
    return [list(path.coords) for path in paths]

# --- Compute total path order ---
def compute_optimal_order(start_pos, path_coords):
    n = len(path_coords)
    indices = list(range(n))
    best_order = None
    best_cost = float('inf')

    for perm in permutations(indices):
        total_cost = 0
        current = start_pos
        for i in perm:
            entry, exit, path = choose_coverage_direction(current, path_coords[i])
            total_cost += euclidean(current, entry)
            total_cost += sum(euclidean(path[j], path[j+1]) for j in range(len(path)-1))
            current = exit
        total_cost += euclidean(current, start_pos)

        if total_cost < best_cost:
            best_cost = total_cost
            best_order = perm

    return best_order

# --- Stitch full path with constrained transitions ---
def stitch_ordered_paths(start_pos, path_coords, order, polygon):
    full_path = []
    current = start_pos

    for i in order:
        entry, exit, path = choose_coverage_direction(current, path_coords[i])
        # Use the original polygon (no buffer) for A* constraint
        transition_path = a_star_path(current, entry, polygon)  # <--- changed here
        full_path.append(transition_path)
        full_path.append(path)
        current = exit

    return_path = a_star_path(current, start_pos, polygon)  # <--- changed here
    full_path.append(return_path)
    return full_path

# --- Visualize full stitched path ---
def visualize_full_path(polygon, convex_cells, full_path_segments):
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    ax.set_title("Drone Full Coverage Path")

    x, y = polygon.exterior.xy
    ax.plot(x, y, 'black', linewidth=2, label='Polygon Boundary')

    for cell in convex_cells:
        x, y = cell.exterior.xy
        ax.plot(x, y, 'grey', linewidth=1, linestyle='--')

    for segment in full_path_segments:
        if len(segment) < 2: continue
        x, y = zip(*segment)
        ax.plot(x, y, color='green', linewidth=1.8)

    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()

def project_to_polygon_boundary(point, polygon):
    """If point is outside polygon, project to nearest point on boundary."""
    pt = Point(point)
    if polygon.contains(pt):
        return point
    # Find nearest point on polygon exterior
    nearest = nearest_points(pt, polygon.exterior)[1]
    return (nearest.x, nearest.y)    

# --- Main block ---
if __name__ == "__main__":
    input_kml_path = r"F:/Meghna/Nidar/buffer/input_polygon.kml"
    output_kml_path = r"optimized_convex_decomposition.kml"
    path_kml_output = "lawnmower_paths.kml"

    spacing = get_path_spacing(
        altitude=12,
        horizontal_fov=130,
        vertical_fov=130,
        overlap=0.3
    )
    spacing_degrees = spacing / 111111
    buffer_space_degrees = 20 / 111111

    polygon = extract_polygon_from_kml(input_kml_path)
    convex_cells = optimize_convex_decomposition(polygon, buffer_space_degrees, spacing_degrees)
    export_convex_cells_to_kml(convex_cells, polygon, output_kml_path)

    paths = generate_lawnmower_path_for_cells(convex_cells, spacing_degrees)
    path_coords = extract_path_coords(paths)

    # User-specified start position (could be outside polygon)
    start_pos = (77.524004, 23.255006)

    # Project start_pos to polygon boundary if outside
    entry_pos = project_to_polygon_boundary(start_pos, polygon)

    # If start_pos is outside, add a transition path from start_pos to entry_pos
    initial_transition = []
    if start_pos != entry_pos:
        initial_transition = [ a_star_path(start_pos, entry_pos, polygon) ]  # <--- changed here

    optimal_order = compute_optimal_order(entry_pos, path_coords)
    full_path_segments = stitch_ordered_paths(entry_pos, path_coords, optimal_order, polygon)

    # Prepend the initial transition if needed
    if initial_transition:
        full_path_segments = initial_transition + full_path_segments

    # --- Save the full stitched path to KML ---
    # Flatten the full_path_segments into LineStrings for KML export
    from shapely.geometry import LineString

    kml_lines = []
    for segment in full_path_segments:
        if len(segment) >= 2:
            kml_lines.append(LineString(segment))

    export_lawnmower_path_to_kml(kml_lines, convex_cells, path_kml_output)

    # Optionally, add a marker for the original start position
    def visualize_full_path_with_start(polygon, convex_cells, full_path_segments, start_pos):
        fig, ax = plt.subplots(1, 1, figsize=(16, 10))
        ax.set_title("Drone Full Coverage Path")

        x, y = polygon.exterior.xy
        ax.plot(x, y, 'black', linewidth=2, label='Polygon Boundary')

        for cell in convex_cells:
            x, y = cell.exterior.xy
            ax.plot(x, y, 'grey', linewidth=1, linestyle='--')

        for segment in full_path_segments:
            if len(segment) < 2: continue
            x, y = zip(*segment)
            ax.plot(x, y, color='green', linewidth=1.8)

        # Mark the original start position
        ax.plot(start_pos[0], start_pos[1], 'ro', markersize=10, label='Original Start')
        # Mark the entry position (projected onto polygon)
        if start_pos != entry_pos:
            ax.plot(entry_pos[0], entry_pos[1], 'bo', markersize=8, label='Entry on Polygon')

        ax.set_aspect('equal')
        ax.grid(True)
        ax.legend()
        plt.tight_layout()
        plt.show()

    visualize_full_path_with_start(polygon, convex_cells, full_path_segments, start_pos)


