"""
Class to read bam files and create a position dictionary
"""

import os
import sys
import pysam
import interlap

class LocusExtractor:
    """
    Class to read bam files and create a position dictionary
    """

    def __init__(self, alignment_file_path, mapping_mode, start_loci_interlaps, stop_loci_interlaps, config, start_output_dict, stop_output_dict):
        self.alignment_file_path = alignment_file_path
        self.mapping_mode = mapping_mode
        self.start_loci_interlaps = start_loci_interlaps
        self.stop_loci_interlaps = stop_loci_interlaps
        self.config = config

        self.start_output_dict = start_output_dict
        self.stop_output_dict = stop_output_dict

        print(
            f"Reading read positions from alignment file: {os.path.basename(alignment_file_path)}"
        )
        self._extract_loci()

    def _init_output_dict(self, output_dict, chrom, read_length):
        """
        Initialize the output dict
        """

        if chrom not in output_dict:
            output_dict[chrom] = {}

        if self.alignment_file_path.stem not in output_dict[chrom]:
            output_dict[chrom][self.alignment_file_path.stem] = {}

        if read_length not in output_dict[chrom][self.alignment_file_path.stem]:
            output_dict[chrom][self.alignment_file_path.stem][read_length] = {}

    def _get_genome_position(self, strand, start, stop):
        """
        Get the genome position
        """

        position_genome = None
        if strand == "-":
            if self.mapping_mode == "threeprime":
                position_genome = start
            elif self.mapping_mode == "fiveprime":
                position_genome = stop
            else:
                sys.exit("Error: Remapping only supported for threeprime or fiveprime mappings.")

        else:
            if self.mapping_mode == "threeprime":
                position_genome = stop
            elif self.mapping_mode == "fiveprime":
                position_genome = start
            else:
                sys.exit("Error: Remapping only supported for threeprime or fiveprime mappings.")
        return position_genome

    def _calculate_window_position(self, strand, position_genome, window, interval):
        """
        Calculate the window position relative to the actual position
        """

        relative_position = position_genome - interval[0]
        if strand == "+":
            position_window = relative_position + window[0]
        else:
            position_window = window[1] - relative_position

        return position_window

    def _fill_output_dict(self, output_dict, position_window, interval, chrom, read_length):
        """
        Fill the output dict
        """

        if position_window not in output_dict[chrom][self.alignment_file_path.stem][read_length]:
            output_dict[chrom][self.alignment_file_path.stem][read_length][position_window] = { interval[2] : 1 }
        else:
            if interval[2] not in output_dict[chrom][self.alignment_file_path.stem][read_length][position_window]:
                output_dict[chrom][self.alignment_file_path.stem][read_length][position_window][interval[2]] = 1

            else:
                output_dict[chrom][self.alignment_file_path.stem][read_length][position_window][interval[2]] += 1

    def _map_reads_to_loci(self, output_dict, loci_interlaps, position_genome, window, chrom, strand, read_length):
        """
        Map the position to the start loci
        """
        if (chrom, strand) not in loci_interlaps:
            return

        for interval in loci_interlaps[(chrom, strand)].find((position_genome, position_genome)):
            position_window = self._calculate_window_position(strand, position_genome, window, interval)

            self._fill_output_dict(output_dict, position_window, interval, chrom, read_length)

    def _extract_loci(self):
        """
        read the alignment file using pysam
        """

        try:
            alignment_file = pysam.AlignmentFile(self.alignment_file_path)
        except (FileNotFoundError, IOError):
            sys.exit("Error: Unable to read alignment file.")
        except ValueError:
            sys.exit("Error: Invalid BAM file format!")

        try:
            for read in alignment_file.fetch():
                chrom = read.reference_name

                if (
                    read.get_tag("NH") > 1
                    or read.mapping_quality < 0
                    or read.is_unmapped
                ):
                    continue

                start = read.reference_start
                read_length = read.query_length # query read length
                stop = start + read_length - 1

                strand = "-" if read.is_reverse else "+"

                self._init_output_dict(self.start_output_dict, chrom, read_length)
                self._init_output_dict(self.stop_output_dict, chrom, read_length)

                position_genome = self._get_genome_position(strand, start, stop)
                self._map_reads_to_loci(self.start_output_dict, self.start_loci_interlaps, position_genome, [-self.config["positionsOutsideORF"], self.config["positionsInsideORF"]], chrom, strand, read_length)
                self._map_reads_to_loci(self.stop_output_dict, self.stop_loci_interlaps, position_genome, [-self.config["positionsInsideORF"], self.config["positionsOutsideORF"]], chrom, strand, read_length)

        except ValueError:
            sys.exit(
                "Error: Ensure that all bam files used \
                 for readcounting have an appropriate index file (.bam.bai). \
                 You can create them using samtools index."
            )


    def output(self):
        """
        Return the output dict.
        """
        return self.start_output_dict, self.stop_output_dict


class IntervalReader():
    """
    Read sam/bam file and create an interlap dictionary for all reads or a selected amount of read lengths
    read_lengths : list (e.g. [30,32,40])
    """
    def __init__(self, alignment_file_path):
        self.alignment_file_path = alignment_file_path

        self.no_accepted_reads_dict = {}
        self.reads_interlap_dict = {}

        print(f"Reading alignment file: {alignment_file_path.stem}")

        self._read_alignment_file()

    def _read_alignment_file(self):
        """
        read the alignment file using pysam
        """

        try:
            alignment_file = pysam.AlignmentFile(self.alignment_file_path)
        except (FileNotFoundError, IOError):
            sys.exit("Error: Unable to read alignment file.")
        except ValueError:
            sys.exit("Error: Invalid BAM file format!")

        tmp_dict = {}
        try:
            for read in alignment_file.fetch():
                chrom = read.reference_name

                if read.get_tag("NH") > 1 or read.mapping_quality < 0 or read.is_unmapped:
                    continue

                start = read.reference_start
                read_length = read.query_length # query read length
                stop = start + read_length - 1

                strand = "-" if read.is_reverse else "+"

                if chrom in self.no_accepted_reads_dict:
                    self.no_accepted_reads_dict[chrom] += 1
                else:
                    self.no_accepted_reads_dict[chrom] = 1

                interval = (start, stop, read_length)

                if (chrom, strand) in tmp_dict:
                    tmp_dict[(chrom, strand)].append(interval)
                else:
                    tmp_dict[(chrom, strand)] = [interval]

            for key, val in tmp_dict.items():
                inter = interlap.InterLap()
                inter.update(val)
                self.reads_interlap_dict[key] = inter

        except ValueError:
            sys.exit("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")

    def output(self):
        return self.reads_interlap_dict, self.no_accepted_reads_dict
