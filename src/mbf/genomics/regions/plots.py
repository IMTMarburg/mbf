# TODO
import hashlib
from typing import List, Optional, Tuple, Union
from pathlib import Path
import dppd
import dppd_plotnine  # noqa: F401
import pypipegraph as ppg
import pandas as pd
import numpy as np
import collections
import mbf.align
import mbf_bam

from mbf.genomics.util import parse_a_or_c
from . import GenomicRegions


dp, X = dppd.dppd()

Sample = mbf.align.lanes._BamDerived

_coverages = {}


def _calc_coverage_bam(gr, sample, flip_anno, flip_column, extend_reads, cache_dir):
    def calc(gr=gr, sample=sample, extend_reads=extend_reads):
        intervals = []
        lens = set()

        for _, row in gr.df.iterrows():
            if flip_column:
                f = row[flip_column]
                if isinstance(f, bool):
                    flip = f
                elif f == "+" or f == 1:
                    flip = False
                elif f == "-" or f == -1:
                    flip = True
                else:
                    raise ValueError("invalid value in flip column", row)
            else:
                flip = False
            lens.add(row["stop"] - row["start"])
            if len(lens) > 1:
                raise ValueError("Different lengths in GR - not supported")
            intervals.append((row["chr"], int(row["start"]), int(row["stop"]), flip))
        bam_filename, index_filename = sample.get_bam_names()
        cov = mbf_bam.calculate_coverage_sum(
            bam_filename, index_filename, intervals, extend_reads
        )
        cov = np.array(cov, dtype=float)
        return cov

    def store(
        coverage, gr_name=gr.name, sample_name=sample.name, extend_reads=extend_reads
    ):
        key = gr_name, sample_name, extend_reads
        if key in _coverages:
            raise ValueError("Douplicate key usage in coverage cache")
        _coverages[key] = coverage

    key = gr.name + "/" + sample.name + "__" + str(extend_reads)
    key = hashlib.sha256(key.encode("utf-8")).hexdigest()
    job = ppg.CachedDataLoadingJob(cache_dir / key, calc, store).depends_on(
        gr.load(), sample.load()
    )
    if flip_anno is not None and flip_column is not None:
        job.depends_on(gr.add_annotator(flip_anno))
    return job


class PlotAveragedCoverage:
    def __init__(
        self,
        gr_to_include: GenomicRegions,
        flip_column: Optional[str],
        extend_reads: int = 0,
    ):
        """
        Plot a (chipseq) signal by averaging the
        signal in @gr_to_include's (equally sized) regions.

        If flip_column is set, the signal in each element is flipped.
        Maybe an annotator

        Names and Colors are a mapping Sample.name => display name | color
        (colors may be color names or #hex. Default is to use scale_many_categories)

        extend_reads: how many bp to extend the read 'downstream'
        """
        if not isinstance(gr_to_include, GenomicRegions):
            raise ValueError("gr_to_include must be a genomic region")
        self.gr_to_include = gr_to_include
        if flip_column is None:
            self.flip_anno = None
            self.flip_column = None
        else:
            self.flip_anno, self.flip_column = parse_a_or_c(flip_column)
        self.extend_reads = int(extend_reads)
        self.samples = []
        self._plot_options = {
            "alpha": 0.5,
            "vertical_bar_pos": "center",
        }
        self.cache_dir = Path("cache/PlotAveragedCoverage/") / gr_to_include.name
        self.cache_dir.mkdir(exist_ok=True, parents=True)
        pass

    def add_sample(
        self,
        sample: Sample,
        name: Optional[str] = None,
        color: Optional[str] = None,
        normalize: Optional[Union[bool, float, Sample]] = None,
    ):
        """Add a sample (e.g. a mbf.align.Lane)
        if @name is None, use sample.name
        if all @color are None, use scale_color_many_categories.
        if only some @color are none: raise

        normalize:
            None | False => no normalization,not even read depth
            True | => Normalize each sample to 1e6 total aligned reads
            Sample =>
                if same genome: substract sample (e.g. IgG)
                else: divide by sample.mapped_reads() / 1e6  (ie. norm by total alignment to other genome)
        """
        if not isinstance(normalize, (bool, Sample, type(None))):
            raise ValueError("invalid value for normalize")
        for x_sample, _, _, x_normalize in self.samples:
            if x_sample == sample and x_normalize == normalize:
                raise ValueError("Same sample with same normalization added again")
        self.samples.append((sample, name, color, normalize))
        self._verify_colors()
        return self

    def _verify_colors(self):
        if not self.samples:
            return True

        seen = set()
        for s, _, c, _ in self.samples:
            seen.add(c)
        if len(seen) == 1 and None in seen:
            return True
        elif None in seen:
            raise ValueError(
                "Mixing defined and None colors - not supported. Either set all, or none"
            )
        elif len(seen) != len(self.samples):
            raise ValueError("Colors duplicated")
        else:
            return True

    def add_samples(
        self,
        sample_tuples: List[
            Tuple[
                Sample,
                Optional[str],
                Optional[str],
                Optional[Union[bool, float, Sample]],
                int,
            ]
        ],
    ):
        """add many samples at once from a list of tuples.
        For tuple contents, see add_sample

        """
        for s, n, c, e in sample_tuples:
            self.add_sample(s, n, c, e)
        return self

    def plot_options(
        self,
        alpha=None,
        vertical_bar_pos=None,
        title=None,
        y_scale_args=None,
        x_scale_args=None,
    ):
        """plotting options.
        vertical_bar_pos may be 'center' (default) to
            have one in the middle (then x axis goes from -bp/2..+bp/2)
        Title may contain %i which get's replaced with the number of regions

        y_scale_args get passed to scale_y_continuous
        """
        if alpha is not None:
            self._plot_options["alpha"] = alpha

        if x_scale_args is not None:
            self._plot_options["x_scale_args"] = x_scale_args
        if y_scale_args is not None:
            self._plot_options["y_scale_args"] = y_scale_args
        if vertical_bar_pos is not None:
            if not isinstance(vertical_bar_pos, int) and vertical_bar_pos not in (
                "center",
                False,
            ):
                raise ValueError(
                    "vertical_bar_pos must be an integer or 'center' or None to disable"
                )
            self._plot_options["vertical_bar_pos"] = vertical_bar_pos
        if title is not None:
            self._plot_options["title"] = title
        return self

    def _calc_coverages(self):
        samples = {x[0].name: x[0] for x in self.samples}
        for sample, _, _, norm in self.samples:
            if isinstance(norm, Sample) and sample.genome == norm.genome:
                samples[norm.name] = norm
        jobs = {}
        for name, sample in samples.items():
            if callable(sample):
                jobs[name] = sample(
                    self.gr_to_include,
                    sample,
                    self.flip_anno,
                    self.flip_column,
                    self.extend_reads,
                    self.cache_dir,
                )

            else:
                jobs[name] = _calc_coverage_bam(
                    self.gr_to_include,
                    sample,
                    self.flip_anno,
                    self.flip_column,
                    self.extend_reads,
                    self.cache_dir,
                )
        return jobs

    def render(self, output_filename):  # noqa: C901
        def calc():
            df = collections.defaultdict(list)
            bp = None
            for sample, _, _, normalization in self.samples:
                cov = _coverages[
                    self.gr_to_include.name, sample.name, self.extend_reads
                ]
                # raise ValueError(sample.name, cov)
                if bp is None:
                    bp = np.arange(cov.shape[0])  # 0...n for the x axis, basically
                if isinstance(normalization, Sample):
                    if normalization.genome == sample.genome:
                        # same genome: *substract* (normalized) reads
                        factor = sample.mapped_reads() / 1e6
                        cov /= factor
                        norm_factor = normalization.mapped_reads() / 1e6
                        cov -= (
                            _coverages[
                                self.gr_to_include.name,
                                normalization.name,
                                self.extend_reads,
                            ]
                            / norm_factor
                        )
                    else:
                        # this assumes you have a spike in of a different species.
                        # the depth in the spike in is obvs very much proportional
                        # to the squencing depth in the sample (they're the same lib!).
                        # so we just use the spike in depth as the normalization depth
                        factor = normalization.mapped_reads() / 1e6
                        cov /= factor
                        # considering the sample depth as well would lead us to
                        # use it twice, essentialy.
                elif isinstance(normalization, float):
                    # normalize by fixed factor
                    cov *= normalization
                elif isinstance(normalization, bool) or normalization is None:
                    if normalization:
                        # true -> normalize to sequening depth
                        factor = sample.mapped_reads() / 1e6
                        cov /= factor
                    else:
                        # use raw read coutns
                        pass
                else:
                    raise NotImplementedError("Unknown normalization")
                # scala is 'per region'
                cov /= len(self.gr_to_include.df)
                df["basepair"].extend(bp)
                df["value"].extend(cov)
                df["sample"].extend([sample.name] * len(bp))
            return pd.DataFrame(df)

        name_repls = {}
        color_map = {}
        names_in_order = []
        for sample, name, color, _norm in self.samples:
            if name:
                name_repls[sample.name] = name
            else:
                name = sample.name
            names_in_order.append(name)
            if color:
                color_map[name] = color

        def plot(df):
            if name_repls:
                df = df.assign(sample=df["sample"].replace(name_repls))

            vertical_bar = self._plot_options["vertical_bar_pos"]
            if vertical_bar == "center":
                df = df.assign(basepair=df["basepair"] - (df["basepair"].max() // 2))

            pdf = dp(df).categorize("sample", names_in_order).p9()

            if vertical_bar == "center":
                pdf = pdf.add_vline(0)
            elif isinstance(vertical_bar, (int, float)):
                pdf = pdf.add_vline(vertical_bar)

            pdf = pdf.add_line(
                "basepair",
                "value",
                color="sample",
                _alpha=self._plot_options["alpha"],
            )
            if color_map:
                pdf = pdf.scale_color_manual(color_map)
            else:
                pdf = pdf.scale_color_many_categories()

            x_scale_args = self._plot_options.get("x_scale_args", None)
            if x_scale_args:
                pdf = pdf.scale_x_continuous(**x_scale_args)

            y_scale_args = self._plot_options.get("y_scale_args", None)
            if y_scale_args:
                pdf = pdf.scale_y_continuous(**y_scale_args)
            title = self._plot_options.get("title", False)
            if title:
                if "%i" in title:
                    title = title % len(self.gr_to_include.df)
                pdf = pdf.title(title)
            # pdf.render(output_filename).pd
            return pdf.pd

        if not hasattr(ppg, "is_ppg2"):
            raise ValueError("Needs pypipegraph2")
        import pypipegraph2 as ppg2

        res = ppg2.PlotJob(output_filename, calc, plot)
        plot_job, (calc_load, calc_cache), _ = res
        calc_cache.depends_on(self._calc_coverages().values())
        calc_cache.depends_on(self.gr_to_include.load())
        plot_job.depends_on_params(
            (
                name_repls,
                color_map,
                names_in_order,
                self._plot_options,
                self.extend_reads,
                sorted(
                    [
                        (x.name, norm.name if hasattr(norm, "name") else str(norm))
                        for x, _, _, norm in self.samples
                    ]
                ),
            )
        )
        if "%i" in self._plot_options.get("title", ""):
            plot_job.depends_on(self.gr_to_include.load())

        return res
