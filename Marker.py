class Marker:
    def __init__(self, well, id, depth = None):
        self.well = well
        self.id = id
        self.correlations = []
        self.depth = depth

    def addCorrel(self, marker):
        """
        :param marker: marker that we want to add to the correlation list of the other marker
        :return:
        """
        id = marker.id
        for m in self.correlations:
            id_m = m.id
            if id == id_m:
                return 0
        self.correlations.append(marker)