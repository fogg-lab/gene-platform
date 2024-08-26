import React from 'react';
import Accordion from '../components/ui/Accordian';
const Guide = () => {
    const faq = [
        {
            question: 'Which databases can be used for the \'Use External Database\' analysis option?',
            answer: 'The GENE Platform supports importing both GEO and GDC datasets.',
        },
        {
            question: 'Why can\'t I find the GEO dataset I\'m looking for?',
            answer: 'Many GEO datasets are not available in a standardized common format, and therefore cannot be parsed without additional configuration.',
        },
        {
            question: 'Why am I not seeing additional the input options that would be necessary for gene set enrichment analysis?',
            answer: 'To perform gene set enrichment analysis, simply navigate to the \'Gene Set Enrichment Analysis\' tab on the Analysis page. By entering into this tab, more input options will become automatically added to the left-hand form.',
        },
        {
            question: 'How do I change the save directory for my analysis runs?',
            answer: 'To change your save directory, navigate to the \'Settings\' icon located on the right side of the navigation bar. There, select \'Change save location\'',
        },
        {
            question: 'What input data file types are accepted for upload?',
            answer: 'The GENE Platform accepts the following file types for counts and coldata:\n\n' +
                    '- Tab-separated values (.tsv or .txt)\n' +
                    '- Comma-separated values (.csv)\n' +
                    '- Gzipped versions of the above (.tsv.gz, .txt.gz, or .csv.gz)',
        },
        {
            question: 'What a cool application!',
            answer: '😎',
        }
    ];
    return (
        <div id="guideContainer">
            <Accordion data={[faq[0]]} />
            <Accordion data={[faq[1]]} />
            <Accordion data={[faq[2]]} />
            <Accordion data={[faq[3]]} />
            <Accordion data={[faq[4]]} />
            <Accordion data={[faq[5]]} />
        </div>
    );
};

export default Guide;
