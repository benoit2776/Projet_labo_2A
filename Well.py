class Well:
    def __init__(self, well_name, nb_markers, markers_depths):
        self.well_name = well_name
        self.markers = []
        self.nb_markers = nb_markers
        self.markers_depths = markers_depths
        self.cost = []
        self.dico_cost = {}

    def get_marker(self, id):
        """
        :param id: id of the marker that we want to get back
        :return: give the marker that corresponds to the id
        """
        for marker in self.markers:
            if marker.id == id:
                return marker
        return None

    def get_well(self, wells_list, well_name):
        """
        :param wells_list: list of wells
        :param well_name: name of the well that we want to get back from the well list
        :return: give the well instance with the name given in argument
        """
        for well in wells_list:
            if well.well_name == well_name:
                return well
        return None

    def marker_presence(self, marker):
        """
        :param marker: marker
        :return: return True if the marker still exists in the marker list
        """
        id = marker.id
        for m in self.markers:
            id_m = m.id
            if id == id_m:
                return True
        return False