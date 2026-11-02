#!/usr/bin/env python3
"""Build supplementary methods tables as a Word document."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT_DIR = REPO_ROOT / "manuscript-outputs" / "tables" / "methods-audit"
OUTPUT_DIR = INPUT_DIR
OUTPUT_DOCX = OUTPUT_DIR / "supplementary_methods_tables.docx"

GRID_CSV = INPUT_DIR / "grid_cell_processing_summary.csv"
SPECIES_CSV = INPUT_DIR / "species_processing_summary.csv"
TRAITS_CSV = INPUT_DIR / "modelled_species_traits_prevalence.csv"

PAGE_WIDTH_IN = 11.0
PAGE_HEIGHT_IN = 8.5
MARGIN_IN = 0.45
TABLE_INDENT_DXA = 120
CONTENT_WIDTH_DXA = int((PAGE_WIDTH_IN - 2 * MARGIN_IN) * 1440)
TABLE_WIDTH_DXA = CONTENT_WIDTH_DXA - TABLE_INDENT_DXA

HEADER_FILL = "E8EEF5"
BORDER_COLOR = "8FA4B8"
HEADING_BLUE = RGBColor(0x2E, 0x74, 0xB5)
HEADING_DARK_BLUE = RGBColor(0x1F, 0x4D, 0x78)
BODY_COLOR = RGBColor(0x00, 0x00, 0x00)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def value(row: dict[str, str], key: str) -> str:
    val = row.get(key, "")
    if val is None or val == "" or val == "NA":
        return ""
    return val


def fmt_int(row: dict[str, str], key: str) -> str:
    raw = value(row, key)
    if raw == "":
        return ""
    return str(int(float(raw)))


def fmt_float(row: dict[str, str], key: str, digits: int = 3) -> str:
    raw = value(row, key)
    if raw == "":
        return ""
    return f"{float(raw):.{digits}f}"


def grid_stage_label(raw: str) -> str:
    labels = {
        "all_grid_cells": "All grid cells",
        "after_land_cover_threshold_pct_lnd_ge_25": "After land-cover threshold",
        "after_effort_filter_good_or_average": "After effort filter",
        "after_effort_filter": "After effort filter",
        "after_missing_covariate_filter_complete_cases": "After missing-covariate filter",
        "after_missing_covariates_filter": "After missing-covariate filter",
        "analysis_grid_cells": "Analysis grid cells",
        "final_model_output_grid_cells": "Final model output grid cells",
    }
    return labels.get(raw, raw.replace("_", " "))


def species_stage_label(raw: str) -> str:
    labels = {
        "recorded_in_atlas_trait_table_occ_gt_0": "Recorded in atlas",
        "removed_by_status_flags": "Removed by status flags",
        "after_status_filter": "After status filter",
        "retained_after_status_flags": "Retained after status flags",
        "recorded_after_grid_covariate_effort_filters_before_min_records": "Recorded after grid, covariate, and effort filters",
        "removed_due_to_insufficient_breeding_records_lt5_iterative": "Removed: insufficient breeding records",
        "removed_insufficient_breeding_records": "Removed: insufficient breeding records",
        "modelled_species_set": "Modelled species set",
        "final_modelled_species": "Final modelled species",
    }
    return labels.get(raw, raw.replace("_", " "))


def set_doc_styles(doc: Document) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = "Calibri"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
    normal.font.size = Pt(11)
    normal.font.color.rgb = BODY_COLOR
    normal.paragraph_format.space_before = Pt(0)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.25

    for style_name, size, color, before, after in [
        ("Heading 1", 16, HEADING_BLUE, 18, 10),
        ("Heading 2", 13, HEADING_BLUE, 14, 7),
        ("Heading 3", 12, HEADING_DARK_BLUE, 10, 5),
    ]:
        style = doc.styles[style_name]
        style.font.name = "Calibri"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
        style.font.size = Pt(size)
        style.font.color.rgb = color
        style.font.bold = True
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.line_spacing = 1.25

    caption = doc.styles.add_style("STableCaption", 1)
    caption.base_style = doc.styles["Normal"]
    caption.font.name = "Calibri"
    caption._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
    caption.font.size = Pt(9.5)
    caption.font.bold = True
    caption.font.color.rgb = HEADING_DARK_BLUE
    caption.paragraph_format.space_before = Pt(4)
    caption.paragraph_format.space_after = Pt(4)
    caption.paragraph_format.line_spacing = 1.1

    note = doc.styles.add_style("TableNote", 1)
    note.base_style = doc.styles["Normal"]
    note.font.name = "Calibri"
    note._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
    note.font.size = Pt(8)
    note.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    note.paragraph_format.space_before = Pt(4)
    note.paragraph_format.space_after = Pt(4)
    note.paragraph_format.line_spacing = 1.1


def configure_page(doc: Document) -> None:
    section = doc.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width = Inches(PAGE_WIDTH_IN)
    section.page_height = Inches(PAGE_HEIGHT_IN)
    section.top_margin = Inches(MARGIN_IN)
    section.bottom_margin = Inches(MARGIN_IN)
    section.left_margin = Inches(MARGIN_IN)
    section.right_margin = Inches(MARGIN_IN)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)


def get_or_add(parent, tag: str):
    existing = parent.find(qn(tag))
    if existing is not None:
        return existing
    element = OxmlElement(tag)
    parent.append(element)
    return element


def remove_children(parent, tag: str) -> None:
    for child in list(parent):
        if child.tag == qn(tag):
            parent.remove(child)


def set_table_geometry(table, widths_in: Iterable[float]) -> None:
    widths_dxa = [int(width * 1440) for width in widths_in]
    table.autofit = False
    table.allow_autofit = False

    tbl = table._tbl
    tbl_pr = tbl.tblPr

    tbl_w = get_or_add(tbl_pr, "w:tblW")
    tbl_w.set(qn("w:w"), str(sum(widths_dxa)))
    tbl_w.set(qn("w:type"), "dxa")

    tbl_ind = get_or_add(tbl_pr, "w:tblInd")
    tbl_ind.set(qn("w:w"), str(TABLE_INDENT_DXA))
    tbl_ind.set(qn("w:type"), "dxa")

    layout = get_or_add(tbl_pr, "w:tblLayout")
    layout.set(qn("w:type"), "fixed")

    cell_mar = get_or_add(tbl_pr, "w:tblCellMar")
    for side, margin in [("top", 80), ("bottom", 80), ("start", 120), ("end", 120)]:
        side_el = get_or_add(cell_mar, f"w:{side}")
        side_el.set(qn("w:w"), str(margin))
        side_el.set(qn("w:type"), "dxa")

    borders = get_or_add(tbl_pr, "w:tblBorders")
    for edge in ["top", "left", "bottom", "right", "insideH", "insideV"]:
        edge_el = get_or_add(borders, f"w:{edge}")
        edge_el.set(qn("w:val"), "single")
        edge_el.set(qn("w:sz"), "4")
        edge_el.set(qn("w:space"), "0")
        edge_el.set(qn("w:color"), BORDER_COLOR)

    tbl_grid = tbl.tblGrid
    if tbl_grid is None:
        tbl_grid = OxmlElement("w:tblGrid")
        tbl.insert(1, tbl_grid)
    remove_children(tbl_grid, "w:gridCol")
    for width in widths_dxa:
        grid_col = OxmlElement("w:gridCol")
        grid_col.set(qn("w:w"), str(width))
        tbl_grid.append(grid_col)

    for row in table.rows:
        for idx, cell in enumerate(row.cells):
            cell.width = Inches(widths_in[idx])
            tc_pr = cell._tc.get_or_add_tcPr()
            tc_w = get_or_add(tc_pr, "w:tcW")
            tc_w.set(qn("w:w"), str(widths_dxa[idx]))
            tc_w.set(qn("w:type"), "dxa")
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP


def set_repeat_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def set_cant_split(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    cant_split = OxmlElement("w:cantSplit")
    tr_pr.append(cant_split)


def shade_cell(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = get_or_add(tc_pr, "w:shd")
    shd.set(qn("w:fill"), fill)


def set_cell_text(cell, text: str, font_size: float, bold: bool = False, align: str = "left") -> None:
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.paragraph_format.line_spacing = 1.0
    if align == "right":
        paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    elif align == "center":
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    else:
        paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = paragraph.add_run(text)
    run.font.name = "Calibri"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Calibri")
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.color.rgb = BODY_COLOR


def add_caption(doc: Document, text: str) -> None:
    paragraph = doc.add_paragraph(style="STableCaption")
    paragraph.add_run(text)


def add_source_note(doc: Document, text: str) -> None:
    paragraph = doc.add_paragraph(style="TableNote")
    paragraph.add_run(text)


def add_table(
    doc: Document,
    rows: list[list[str]],
    headers: list[str],
    widths_in: list[float],
    numeric_columns: set[int] | None = None,
    body_font: float = 7.0,
    header_font: float = 7.0,
) -> None:
    numeric_columns = numeric_columns or set()
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    set_table_geometry(table, widths_in)
    set_repeat_header(table.rows[0])
    set_cant_split(table.rows[0])

    for idx, header in enumerate(headers):
        cell = table.rows[0].cells[idx]
        shade_cell(cell, HEADER_FILL)
        align = "right" if idx in numeric_columns else "left"
        set_cell_text(cell, header, header_font, bold=True, align=align)

    for row_values in rows:
        row = table.add_row()
        set_cant_split(row)
        for idx, text in enumerate(row_values):
            cell = row.cells[idx]
            align = "right" if idx in numeric_columns else "left"
            set_cell_text(cell, text, body_font, align=align)

    set_table_geometry(table, widths_in)


def build_grid_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    output = []
    for row in rows:
        output.append(
            [
                fmt_int(row, "atlas"),
                value(row, "period"),
                fmt_int(row, "total_starting_grid_cells"),
                fmt_int(row, "filtered_area"),
                fmt_int(row, "filtered_incomplete_environmental_variables"),
                fmt_int(row, "filtered_effort"),
                fmt_int(row, "grid_cells_retained"),
            ]
        )
    return output


def build_species_processing_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    output = []
    for row in rows:
        output.append(
            [
                fmt_int(row, "total_species_across_all_atlases_before_filtering"),
                fmt_int(row, "filtered_uncertain_breeding"),
                fmt_int(row, "filtered_introduced_domestic"),
                fmt_int(row, "filtered_rarity_min_records"),
                fmt_int(row, "species_retained"),
            ]
        )
    return output


def build_trait_rows(rows: list[dict[str, str]]) -> list[list[str]]:
    sorted_rows = sorted(rows, key=lambda row: value(row, "latin_name_dof").lower())
    output = []
    for row in sorted_rows:
        output.append(
            [
                value(row, "latin_name_dof"),
                value(row, "english_DOF"),
                fmt_float(row, "species_thermal_index", digits=2),
                value(row, "Migration_a3_DOF"),
                value(row, "foraging_guild_consensus"),
                fmt_float(row, "atlas1_prevalence", digits=3),
                fmt_float(row, "atlas2_prevalence", digits=3),
                fmt_float(row, "atlas3_prevalence", digits=3),
            ]
        )
    return output


def main() -> None:
    grid_rows = read_csv(GRID_CSV)
    species_rows = read_csv(SPECIES_CSV)
    trait_rows = read_csv(TRAITS_CSV)

    doc = Document()
    configure_page(doc)
    set_doc_styles(doc)

    title = doc.add_paragraph(style="Heading 1")
    title.add_run("Supplementary Methods Tables")
    add_source_note(
        doc,
        "Source CSVs: grid_cell_processing_summary.csv, species_processing_summary.csv, and modelled_species_traits_prevalence.csv in manuscript-outputs/tables/methods-audit/.",
    )

    add_caption(
        doc,
        "S Table 1. Grid cell processing by atlas period.",
    )
    add_table(
        doc,
        build_grid_rows(grid_rows),
        [
            "Atlas",
            "Period",
            "Total starting grid cells",
            "# filtered (area)",
            "# filtered (incomplete environmental variables)",
            "# filtered (effort)",
            "# grid cells retained",
        ],
        [0.55, 0.80, 1.65, 1.25, 2.35, 1.25, 1.32],
        numeric_columns={0, 2, 3, 4, 5, 6},
        body_font=8.0,
        header_font=7.6,
    )
    add_source_note(
        doc,
        "Generated from grid_cell_processing_summary.csv. Area filtering is the >=25% land-area threshold; retained cells are the unique fitted-model site levels.",
    )

    doc.add_page_break()
    add_caption(
        doc,
        "S Table 2. Species processing across all atlas periods.",
    )
    add_table(
        doc,
        build_species_processing_rows(species_rows),
        [
            "Total species across all atlases (before filtering)",
            "# filtered (uncertain breeding)",
            "# filtered (introduced/domestic)",
            "# filtered (rarity; <5 records in at least one atlas)",
            "# species retained",
        ],
        [2.50, 1.80, 1.80, 2.60, 1.42],
        numeric_columns={0, 1, 2, 3, 4},
        body_font=8.2,
        header_font=7.8,
    )
    add_source_note(
        doc,
        "Generated from species_processing_summary.csv. The rarity column combines explicit rare-species flags with the final minimum-record threshold used for the modelled species set.",
    )

    doc.add_page_break()
    add_caption(
        doc,
        "S Table 3. Modelled species traits and atlas prevalence.",
    )
    add_table(
        doc,
        build_trait_rows(trait_rows),
        [
            "Latin name",
            "English name",
            "STI",
            "Migration",
            "Foraging guild",
            "A1 prev.",
            "A2 prev.",
            "A3 prev.",
        ],
        [1.62, 1.62, 0.52, 1.35, 2.06, 0.92, 0.92, 0.92],
        numeric_columns={2, 5, 6, 7},
        body_font=6.5,
        header_font=6.5,
    )
    add_source_note(
        doc,
        "Generated from modelled_species_traits_prevalence.csv. Latin names follow DOF names, including subspecies names where present; STI = species thermal index; prev. = breeding-record prevalence among retained analysis grid cells.",
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    doc.save(OUTPUT_DOCX)

    print(f"wrote={OUTPUT_DOCX}")
    print(f"grid_rows={len(grid_rows)}")
    print(f"species_processing_rows={len(species_rows)}")
    print(f"species_trait_rows={len(trait_rows)}")


if __name__ == "__main__":
    main()
