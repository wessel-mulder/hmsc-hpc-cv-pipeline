from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


PROJECT_ROOT = Path(__file__).resolve().parents[3]
FULL_VS_CV_FIGURE_PATH = (
    PROJECT_ROOT
    / "manuscript-outputs"
    / "figures"
    / "supplement"
    / "2026-07-27-sfig-full-vs-cv-accuracy.png"
)
PREVALENCE_FIGURE_PATH = (
    PROJECT_ROOT
    / "manuscript-outputs"
    / "figures"
    / "supplement"
    / "2026-07-27-sfig2-prevalence-vs-cv-auc-tjur.png"
)
COMPONENT_INFLUENCE_FIGURE_PATH = (
    PROJECT_ROOT
    / "manuscript-outputs"
    / "figures"
    / "supplement"
    / "2026-07-27-sfig-component-influence-sankey-bump.png"
)
OUTPUT_PATH = (
    PROJECT_ROOT
    / "manuscript-outputs"
    / "supplementary-material"
    / "supplementary-figure-1-full-vs-cv-accuracy.docx"
)


def set_font(run, name="Calibri"):
    run.font.name = name
    run._element.rPr.rFonts.set(qn("w:ascii"), name)
    run._element.rPr.rFonts.set(qn("w:hAnsi"), name)


def set_paragraph_spacing(paragraph, before=0, after=6, line_spacing=1.10):
    paragraph.paragraph_format.space_before = Pt(before)
    paragraph.paragraph_format.space_after = Pt(after)
    paragraph.paragraph_format.line_spacing = line_spacing


def apply_document_styles(doc):
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1)
    section.right_margin = Inches(1)
    section.bottom_margin = Inches(1)
    section.left_margin = Inches(1)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)

    normal = doc.styles["Normal"]
    normal.font.name = "Calibri"
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
    normal.font.size = Pt(11)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    normal.paragraph_format.space_before = Pt(0)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.10

    heading_1 = doc.styles["Heading 1"]
    heading_1.font.name = "Calibri"
    heading_1._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
    heading_1._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
    heading_1.font.size = Pt(16)
    heading_1.font.bold = True
    heading_1.font.color.rgb = RGBColor(46, 116, 181)
    heading_1.paragraph_format.space_before = Pt(16)
    heading_1.paragraph_format.space_after = Pt(8)
    heading_1.paragraph_format.line_spacing = 1.10


def set_picture_alt_text(paragraph, title, description):
    drawing_nodes = paragraph._p.xpath(".//wp:docPr")
    for drawing_node in drawing_nodes:
        drawing_node.set("title", title)
        drawing_node.set("descr", description)


def add_supplementary_figure(
    doc,
    figure_number,
    figure_path,
    heading,
    alt_title,
    alt_description,
    caption_body,
):
    title = doc.add_paragraph(style="Heading 1")
    title.alignment = WD_ALIGN_PARAGRAPH.LEFT
    title.add_run(heading)

    figure = doc.add_paragraph()
    figure.alignment = WD_ALIGN_PARAGRAPH.CENTER
    set_paragraph_spacing(figure, before=4, after=8, line_spacing=1.0)
    figure_run = figure.add_run()
    figure_run.add_picture(str(figure_path), width=Inches(6.5))
    set_picture_alt_text(
        figure,
        alt_title,
        alt_description,
    )

    caption = doc.add_paragraph()
    set_paragraph_spacing(caption, before=4, after=6, line_spacing=1.10)
    caption_run = caption.add_run(f"Supplementary Figure {figure_number}. ")
    caption_run.bold = True
    set_font(caption_run)
    body_run = caption.add_run(caption_body)
    set_font(body_run)


def build_docx():
    missing_figures = [
        figure_path
        for figure_path in (
            FULL_VS_CV_FIGURE_PATH,
            PREVALENCE_FIGURE_PATH,
            COMPONENT_INFLUENCE_FIGURE_PATH,
        )
        if not figure_path.exists()
    ]
    if missing_figures:
        raise FileNotFoundError(
            "Missing figure PNG(s): " + ", ".join(str(path) for path in missing_figures)
        )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    doc = Document()
    apply_document_styles(doc)

    add_supplementary_figure(
        doc=doc,
        figure_number=1,
        figure_path=FULL_VS_CV_FIGURE_PATH,
        heading="Supplementary Figure 1",
        alt_title="Supplementary Figure 1. Full-model and cross-validated accuracy",
        alt_description=(
            "Violin, point, and boxplot panels comparing full-model and "
            "10-fold cross-validated AUC, Tjur R2, and RMSE across the "
            "1970s, 1990s, and 2010s atlas periods."
        ),
        caption_body=(
            "Full-model and 10-fold cross-validated model accuracy across atlas periods. "
            "Each panel shows species-level accuracy values for AUC, Tjur R2, and RMSE; "
            "diamonds and text labels indicate median values."
        ),
    )

    doc.add_page_break()

    add_supplementary_figure(
        doc=doc,
        figure_number=2,
        figure_path=PREVALENCE_FIGURE_PATH,
        heading="Supplementary Figure 2",
        alt_title=(
            "Supplementary Figure 2. Species prevalence and cross-validated "
            "AUC and Tjur R2"
        ),
        alt_description=(
            "Scatterplots showing species prevalence against 10-fold "
            "cross-validated AUC and Tjur R2. Points are species-period "
            "combinations coloured by atlas period, and black curves are "
            "quadratic linear-model fits with R2 labels."
        ),
        caption_body=(
            "Species prevalence and 10-fold cross-validated AUC and Tjur R2. "
            "Points are species-period combinations coloured by atlas period. "
            "Black curves show quadratic linear-model fits, with facet-specific "
            "R2 values printed in each panel."
        ),
    )

    doc.add_page_break()

    add_supplementary_figure(
        doc=doc,
        figure_number=3,
        figure_path=COMPONENT_INFLUENCE_FIGURE_PATH,
        heading="Supplementary Figure 3",
        alt_title=(
            "Supplementary Figure 3. Component influence through time"
        ),
        alt_description=(
            "Alluvial sankey-bump plot showing the mean absolute variance "
            "explained by the spatial, land-use, and climate components in "
            "the 1970s, 1990s, and 2010s atlas periods."
        ),
        caption_body=(
            "Temporal change in the mean absolute variance explained by three "
            "aggregated components. The spatial component is the site-level "
            "random effect; the land-use component combines land-use "
            "heterogeneity and urban, cropland, pasture, forest, and "
            "grass/shrubland coverage; and the climate component combines "
            "breeding-season temperature and precipitation. Values are based "
            "on Tjur R2-scaled variance partitioning across the 157 modelled "
            "species."
        ),
    )

    doc.save(OUTPUT_PATH)
    print(OUTPUT_PATH)


if __name__ == "__main__":
    build_docx()
