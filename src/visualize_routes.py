#!/usr/bin/env python3
"""
Visualize saved aizynthfinder output JSON as route images.

Usage:
    conda run -n retro python src/visualize_routes.py \
        --input outputs/26-06-12-hesperidin_combined_curated.json \
        --outdir outputs/route_images
"""
import argparse
import json
import os
from pathlib import Path

from aizynthfinder.reactiontree import ReactionTree


def visualize(input_file: str, outdir: str, fmt: str = "png") -> None:
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    with open(input_file) as f:
        routes = json.load(f)

    if not isinstance(routes, list):
        raise ValueError(f"Expected JSON list of routes, got {type(routes)}")

    for i, route_dict in enumerate(routes):
        rt = ReactionTree.from_dict(route_dict)
        img = rt.to_image()
        nsteps = len(list(rt.reactions()))
        outfile = outpath / f"route_{i + 1:02d}.{fmt}"
        img.save(outfile)
        print(
            f"Route {i + 1}: solved={rt.is_solved}, steps={nsteps}, "
            f"size={img.size} -> {outfile}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Render aizynthfinder output JSON routes to images"
    )
    parser.add_argument("--input", required=True, help="Path to output JSON")
    parser.add_argument("--outdir", default="outputs/route_images", help="Output dir")
    parser.add_argument("--fmt", default="png", choices=["png", "svg"], help="Format")
    args = parser.parse_args()
    visualize(args.input, args.outdir, args.fmt)
