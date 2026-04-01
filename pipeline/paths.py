from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DEFAULT_CONTINENTS = ("af", "eu", "na", "sa", "as", "oc")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_results_root(path_like: str | os.PathLike[str]) -> Path:
    path = Path(path_like).expanduser()
    if not path.is_absolute():
        path = (repo_root() / path).resolve()
    if path.name == "results":
        return path
    return path / "results"


def _resolve_path(value: str | None, *, base_dir: Path) -> Path | None:
    if value in (None, ""):
        return None

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


@dataclass(frozen=True)
class ProjectPaths:
    repo_root: Path
    config_path: Path | None
    data_root: Path | None
    results_root: Path
    swot_vector_dir: Path
    swot_nodes_dir: Path
    input_created_root: Path
    fabdem_vrt: Path

    @property
    def new_segments_root(self) -> Path:
        return self.results_root / "new_segments"

    @property
    def new_segments_vector_dir(self) -> Path:
        return self.new_segments_root / "vector"

    @property
    def new_segments_node_dir(self) -> Path:
        return self.new_segments_root / "node"

    @property
    def new_segments_vector_cont_dir(self) -> Path:
        return self.new_segments_root / "vector_cont"

    @property
    def reference_tables_dir(self) -> Path:
        return self.results_root / "reference_tables"

    @property
    def all_dir(self) -> Path:
        return self.results_root / "all"

    @property
    def single_values_dir(self) -> Path:
        return self.results_root / "single_values"

    @property
    def reach_averaged_dir(self) -> Path:
        return self.results_root / "reach_averaged"

    @property
    def single_smoothed_dir(self) -> Path:
        return self.results_root / "single_smoothed"

    @property
    def centerline_dir(self) -> Path:
        return self.results_root / "centerline"

    @property
    def cycles_dir(self) -> Path:
        return self.results_root / "cycles"

    def ensure_step1_dirs(self) -> None:
        for path in (
            self.results_root,
            self.new_segments_root,
            self.new_segments_vector_dir,
            self.new_segments_node_dir,
            self.new_segments_vector_cont_dir,
            self.reference_tables_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def ensure_step2_dirs(self) -> None:
        for path in (
            self.results_root,
            self.all_dir,
            self.centerline_dir,
            self.cycles_dir,
            self.input_created_root,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def validate_step1_inputs(self) -> None:
        missing = []
        if not self.swot_vector_dir.exists():
            missing.append(f"swot_vector_dir: {self.swot_vector_dir}")
        if not self.swot_nodes_dir.exists():
            missing.append(f"swot_nodes_dir: {self.swot_nodes_dir}")

        if missing:
            raise FileNotFoundError(
                "Step 1 inputs are missing. Configure the paths in "
                f"{self.repo_root / 'config' / 'paths.local.json'} or create the expected directories:\n"
                + "\n".join(missing)
            )


def load_project_paths(config_path: str | os.PathLike[str] | None = None) -> ProjectPaths:
    root = repo_root()
    config_dir = root / "config"

    selected_config = None
    if config_path is not None:
        selected_config = Path(config_path).expanduser()
        if not selected_config.is_absolute():
            selected_config = (root / selected_config).resolve()
    else:
        env_config = os.getenv("RIVER_CONFINEMENT_PATHS")
        if env_config:
            selected_config = Path(env_config).expanduser()
            if not selected_config.is_absolute():
                selected_config = (root / selected_config).resolve()
        else:
            local_config = config_dir / "paths.local.json"
            if local_config.exists():
                selected_config = local_config

    raw: dict[str, Any] = {}
    if selected_config is not None:
        if not selected_config.exists():
            raise FileNotFoundError(f"Config file not found: {selected_config}")
        raw = _load_json(selected_config)

    config_base = selected_config.parent if selected_config is not None else root

    data_root = _resolve_path(raw.get("data_root"), base_dir=config_base)
    if data_root is None:
        local_data_root = root / "data"
        if local_data_root.exists():
            data_root = local_data_root

    results_root = _resolve_path(raw.get("results_root"), base_dir=config_base)
    if results_root is None:
        results_root = root / "results"

    input_created_root = _resolve_path(raw.get("input_created_root"), base_dir=config_base)
    if input_created_root is None:
        input_created_root = root / "input_created"

    swot_vector_dir = _resolve_path(raw.get("swot_vector_dir"), base_dir=config_base)
    if swot_vector_dir is None and data_root is not None:
        swot_vector_dir = data_root / "SWOT_vector"
    if swot_vector_dir is None:
        swot_vector_dir = root / "input" / "SWOT_vector"

    swot_nodes_dir = _resolve_path(raw.get("swot_nodes_dir"), base_dir=config_base)
    if swot_nodes_dir is None and data_root is not None:
        swot_nodes_dir = data_root / "SWOT_nodes"
    if swot_nodes_dir is None:
        swot_nodes_dir = root / "input" / "SWOT_nodes"

    fabdem_vrt = _resolve_path(raw.get("fabdem_vrt"), base_dir=config_base)
    if fabdem_vrt is None:
        fabdem_vrt = input_created_root / "FAB_dem_vrt.vrt"

    return ProjectPaths(
        repo_root=root,
        config_path=selected_config,
        data_root=data_root,
        results_root=results_root.resolve(),
        swot_vector_dir=swot_vector_dir.resolve(),
        swot_nodes_dir=swot_nodes_dir.resolve(),
        input_created_root=input_created_root.resolve(),
        fabdem_vrt=fabdem_vrt.resolve(),
    )
