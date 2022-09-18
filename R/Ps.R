#' Compute a distance-dependent contact frequency, a.k.a. P(s)
#' 
#' P(s) will be approximated if no pairs are provided, or the exact P(s) 
#' will be computed if a `.pairs` file is added to the `contacts` object 
#' using `pairsFile(x) <- "..."`. 
#' 
#' @rdname Ps
#' 
#' @param x A `contacts` object
#' @param by_chr by_chr
#' @param filtered_chr filtered_chr
#' @return a tibble
#'
#' @import tibble
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr group_split
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select
#' @importFrom dplyr tally
#' @importFrom dplyr arrange
#' @importFrom dplyr lead
#' @export
#' @examples 
#' library(HiContacts)
#' data(contacts_yeast)
#' getPs(contacts_yeast)

getPs <- function(
    x, 
    by_chr = FALSE, 
    filtered_chr = c('XII', 'chrXII', 'chr12', '12', 'Mito', 'MT', 'chrM')
) {
    pairsFile <- pairsFile(x)
    if (is.null(pairsFile)) {
        # stop("Please provide a pairsFile for `x`. Aborting now.")
        message("pairsFile not specified. The P(s) curve will be an approximation.")
        pairs <- scores(x, 'raw')
        df <- tibble::tibble(
            chr = as.vector(GenomeInfoDb::seqnames(InteractionSet::anchors(pairs)[[1]])),
            distance = InteractionSet::pairdist(pairs, type = 'gap'),
            n = pairs$score
        ) %>% 
            tidyr::drop_na() %>% 
            dplyr::filter(!chr %in% filtered_chr) %>% 
            dplyr::mutate(binned_distance = PsBreaks()$break_pos[findInterval(distance, vec = PsBreaks()$break_pos, all.inside = TRUE)])
        if (by_chr) {
            df <- dplyr::group_by(df, chr, binned_distance)
        } 
        else {
            df <- dplyr::group_by(df, binned_distance)
        }
        d <- dplyr::summarize(df, ninter = sum(n)) %>%
            dplyr::mutate(p = ninter/sum(ninter)) %>% 
            dplyr::left_join(PsBreaks(), by = c('binned_distance' = 'break_pos')) %>% 
            dplyr::mutate(norm_p = p / binwidth)
        if (by_chr) {
            d <- dplyr::group_by(d, chr)
        } 
        else {
            d <- d
        }
        ps <- dplyr::group_split(d) %>% 
            lapply(function(x) {
                x %>% 
                    dplyr::mutate(
                        norm_p_unity = norm_p / 
                            {dplyr::slice(x, which.min(abs(x$binned_distance - 100000))) %>% dplyr::pull(norm_p)}
                    ) %>% 
                    dplyr::mutate(
                        slope = (log10(dplyr::lead(norm_p)) - log10(norm_p)) / 
                            (log10(dplyr::lead(binned_distance)) - log10(binned_distance))
                    ) %>% 
                    dplyr::mutate(
                        slope = c(0, predict(
                            loess(slope ~ binned_distance, span = 0.5, data = .)
                        ))
                    )
            }) %>% 
            dplyr::bind_rows()
        if (by_chr) {
            ps <- dplyr::select(ps, chr, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(chr, binned_distance)
        } 
        else {
            ps <- dplyr::select(ps, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(binned_distance)
        }
        return(ps)
    }
    else {
        message("Importing pairs file ", pairsFile, " in memory. This may take a while...")
        pairs <- pairs2gi(pairsFile)
        df <- tibble::tibble(
            chr = as.vector(GenomeInfoDb::seqnames(InteractionSet::anchors(pairs)[[1]])),
            distance = pairs$distance
        ) %>% 
            tidyr::drop_na() %>% 
            dplyr::filter(!chr %in% filtered_chr) %>% 
            dplyr::mutate(binned_distance = PsBreaks()$break_pos[findInterval(distance, vec = PsBreaks()$break_pos, all.inside = TRUE)])
        if (by_chr) {
            df <- dplyr::group_by(df, chr, binned_distance)
        } 
        else {
            df <- dplyr::group_by(df, binned_distance)
        }
        d <- dplyr::tally(df, name = 'ninter') %>%
            dplyr::mutate(p = ninter/sum(ninter)) %>% 
            dplyr::left_join(PsBreaks(), by = c('binned_distance' = 'break_pos')) %>% 
            dplyr::mutate(norm_p = p / binwidth)
        if (by_chr) {
            d <- dplyr::group_by(d, chr)
        } 
        else {
            d <- d
        }
        ps <- dplyr::group_split(d) %>% 
            lapply(function(x) {
                dplyr::mutate(x, norm_p_unity = norm_p / {dplyr::slice(x, which.min(abs(x$binned_distance - 100000))) %>% dplyr::pull(norm_p)}) %>% 
                dplyr::mutate(slope = (log10(dplyr::lead(norm_p)) - log10(norm_p)) / (log10(dplyr::lead(binned_distance)) - log10(binned_distance))) %>% 
                dplyr::mutate(slope = c(0, predict(loess(slope ~ binned_distance, span = 0.5, data = .))))
            }) %>% 
            dplyr::bind_rows()
        if (by_chr) {
            ps <- dplyr::select(ps, chr, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(binned_distance)
        } 
        else {
            ps <- dplyr::select(ps, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(binned_distance)
        }
        return(ps)
    }
}

#' @rdname Ps
#' 
#' @param ... ...
#' @param xlim xlim
#' @param ylim ylim
#' @return ggplot
#'
#' @import ggplot2
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @export
#' @examples 
#' ## Single P(s)
#' 
#' library(HiContacts)
#' data(contacts_yeast)
#' ps <- getPs(contacts_yeast)
#' plotPs(ps, ggplot2::aes(x = binned_distance, y = norm_p))
#' 
#' ## Comparing several P(s)
#' 
#' library(HiContacts)
#' data(contacts_yeast)
#' data(contacts_yeast_eco1)
#' ps_wt <- getPs(contacts_yeast)
#' ps_wt$sample <- 'WT'
#' ps_eco1 <- getPs(contacts_yeast_eco1)
#' ps_eco1$sample <- 'eco1'
#' ps <- rbind(ps_wt, ps_eco1)
#' plotPs(ps, ggplot2::aes(x = binned_distance, y = norm_p, group = sample, color = sample))

plotPs <- function(..., xlim = c(5000, 4.99e5), ylim = c(1e-8, 1e-4)) {
    gg <- ggplot2::ggplot(...) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_y_log10(
            limits = ylim, 
            expand = c(0, 0), 
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        ggplot2::scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        ) + 
        ggplot2::annotation_logticks() + 
        ggplot2::labs(x = "Genomic distance", y = "Contact frequency")
    gg
}

#' @rdname Ps
#' 
#' @return ggplot
#' 
#' @import ggplot2
#' @export
#' @examples 
#' library(HiContacts)
#' data(contacts_yeast)
#' ps <- getPs(contacts_yeast)
#' plotPsSlope(ps, ggplot2::aes(x = binned_distance, y = slope))

plotPsSlope <- function(..., xlim = c(5000, 4.99e5), ylim = c(-3, 0)) {
    gg <- ggplot2::ggplot(...) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_y_continuous(
            limits = ylim, 
            expand = c(0, 0)
        ) +
        ggplot2::scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        ) + 
        ggplot2::annotation_logticks(sides = "b") + 
        ggplot2::labs(x = "Genomic distance", y = "Slope of P(s)")
    gg
}

#' @rdname Ps
#'
#' @return tbl

PsBreaks <- function() {
    structure(list(break_pos = c(
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 16, 17, 19, 21, 23, 26, 28, 31, 34, 37, 41, 45, 50,
        55, 60, 66, 73, 80, 88, 97, 107, 117, 129, 142, 156, 172, 189,
        208, 229, 252, 277, 304, 335, 368, 405, 446, 490, 539, 593, 653,
        718, 790, 869, 956, 1051, 1156, 1272, 1399, 1539, 1693, 1862,
        2048, 2253, 2479, 2726, 2999, 3299, 3629, 3992, 4391, 4830, 5313,
        5844, 6429, 7072, 7779, 8557, 9412, 10354, 11389, 12528, 13781,
        15159, 16675, 18342, 20176, 22194, 24413, 26855, 29540, 32494,
        35743, 39318, 43249, 47574, 52332, 57565, 63322, 69654, 76619,
        84281, 92709, 101980, 112178, 123396, 135735, 149309, 164240,
        180664, 198730, 218603, 240463, 264510, 290961, 320057, 352063,
        387269, 425996, 468595, 515455, 567000, 623700, 686070, 754677,
        830145, 913160, 1004475, 1104923, 1215415, 1336957, 1470653,
        1617718, 1779490, 1957439, 2153182, 2368501, 2605351, 2865886,
        3152474, 3467722, 3814494, 4195943, 4615538, 5077092, 5584801,
        6143281, 6757609, 7433370, 8176707, 8994377, 9893815, 10883197,
        11971516, 13168668, 14485535, 15934088, 17527497, 19280247, 21208271,
        23329099, 25662008, 28228209, 31051030, 34156133, 37571747, 41328921,
        45461813, 50007995, 55008794, 60509674, 66560641, 73216705, 80538375,
        88592213, 97451434, 107196578, 117916236, 129707859, 142678645,
        156946509, 172641160, 189905276, 208895804, 229785385, 252763923,
        278040315, 305844347, 336428781, 370071660, 407078826, 447786708,
        492565379, 541821917, 596004109, 655604519, 721164971, 793281468,
        872609615, 959870577, 1055857635, 1161443398
    ), binwidth = c(
        1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 3, 2, 3, 3,
        3, 4, 4, 5, 5, 5, 6, 7, 7, 8, 9, 10, 10, 12, 13, 14, 16, 17,
        19, 21, 23, 25, 27, 31, 33, 37, 41, 44, 49, 54, 60, 65, 72, 79,
        87, 95, 105, 116, 127, 140, 154, 169, 186, 205, 226, 247, 273,
        300, 330, 363, 399, 439, 483, 531, 585, 643, 707, 778, 855, 942,
        1035, 1139, 1253, 1378, 1516, 1667, 1834, 2018, 2219, 2442, 2685,
        2954, 3249, 3575, 3931, 4325, 4758, 5233, 5757, 6332, 6965, 7662,
        8428, 9271, 10198, 11218, 12339, 13574, 14931, 16424, 18066,
        19873, 21860, 24047, 26451, 29096, 32006, 35206, 38727, 42599,
        46860, 51545, 56700, 62370, 68607, 75468, 83015, 91315, 100448,
        110492, 121542, 133696, 147065, 161772, 177949, 195743, 215319,
        236850, 260535, 286588, 315248, 346772, 381449, 419595, 461554,
        507709, 558480, 614328, 675761, 743337, 817670, 899438, 989382,
        1088319, 1197152, 1316867, 1448553, 1593409, 1752750, 1928024,
        2120828, 2332909, 2566201, 2822821, 3105103, 3415614, 3757174,
        4132892, 4546182, 5000799, 5500880, 6050967, 6656064, 7321670,
        8053838, 8859221, 9745144, 10719658, 11791623, 12970786, 14267864,
        15694651, 17264116, 18990528, 20889581, 22978538, 25276392, 27804032,
        30584434, 33642879, 37007166, 40707882, 44778671, 49256538, 54182192,
        59600410, 65560452, 72116497, 79328147, 87260962, 95987058, 105585763,
        116144340
    )), row.names = c(NA, -205L), class = c(
        "tbl_df", "tbl",
        "data.frame"
    ))
}