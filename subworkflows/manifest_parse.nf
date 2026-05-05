//
// Check input samplesheet and get read channels
//

workflow MANIFEST_PARSE {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel
        .fromPath( samplesheet )
        .ifEmpty {error("Cannot find path file ${samplesheet}")}
        .splitCsv ( header:true, sep:',' )
        .map { create_assembly_channels(it) }
        .set { assemblies }

    emit:
    assemblies
    // pre_generated_annotation_channel for future
}

// Function to get list of [ meta, assembly ]
def create_assembly_channels(LinkedHashMap row) {
    def array = []
    def meta = [:]

    // validate required columns
    if ( !row.containsKey('ID') ) {
        error("ERROR: Manifest is missing required column 'ID'")
    }
    if ( !row.containsKey('assembly') ) {
        error("ERROR: Manifest is missing required column 'assembly'")
    }

    //for bakta
    meta.ID = row.ID.replace('#', '_')

    // check short reads
    if ( !file(row.assembly).exists() ) {
        error("ERROR: Please check input samplesheet -> Assembly file does not exist!\n${row.assembly}")
    }

    def assembly = file(row.assembly)

    array = [ meta, assembly ]
    return array
}

/* Function to get list of [ meta, assembly ]
def create_annotations_channels(LinkedHashMap row) {
    def meta = [:]

    //for bakta
    meta.ID = row.ID.replace('#', '_')

    def array = []
    
    // check short reads
    if ( row.annotations ) {
        annotations = file(row.annotations)
    }

    array = [ meta, annotations ]
    return array
}
*/