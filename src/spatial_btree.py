import pylinear.array as num
import pytools




def do_boxes_intersect((bl1,tr1), (bl2,tr2)):
    (dimension,) = bl1.shape
    for i in range(0, dimension):
        if max(bl1[i], bl2[i]) > min(tr1[i], tr2[i]):
            return False
    return True




def _get_elements_bounding_box(elements):
    if len(elements) == 0:
        raise RuntimeError, "Cannot get the bounding box of no elements."

    bboxes = [ box for el,box in elements ]
    bottom_lefts = [ bl for bl,tr in bboxes ]
    top_rights = [ tr for bl,tr in bboxes ]
    return reduce(num.minimum, bottom_lefts), reduce(num.maximum, top_rights)





def make_buckets(bottom_left, top_right, allbuckets):
    (dimensions,) = bottom_left.shape

    half = (top_right - bottom_left) / 2.
    def do(dimension, pos):
        if dimension == dimensions:
            origin = bottom_left + num.multiply(pos, half)
            bucket = SpatialBinaryTreeBucket(origin, origin + half)
            allbuckets.append(bucket)
            return bucket
        else:
            pos[ dimension ] = 0
            first = do(dimension + 1, pos)
            pos[ dimension ] = 1
            second = do(dimension + 1, pos)
            return [ first, second ]

    return do(0, num.zeros((dimensions,), num.Float))




class SpatialBinaryTreeBucket:
    """This class represents one bucket in a spatial binary tree.
    It automatically decides whether it needs to create more subdivisions 
    beneath itself or not.
    """

    def __init__(self, bottom_left, top_right):
        self.BottomLeft = bottom_left
        self.TopRight = top_right
        self.Center = (bottom_left + top_right) / 2
        
        # As long as Buckets is None, there are no subdivisions
        self.Buckets = None
        self.Elements = []
    
    def bottom_left(self):
        return self.BottomLeft

    def top_right(self):
        return self.TopRight

    def add_element(self, element, bbox):
        (dimensions,) = self.BottomLeft.shape
        if self.Buckets is None:
            # No subdivisions yet.
            if len(self.Elements) > 8 * 2**dimensions:
                # Too many elements. Need to subdivide.
                self.AllBuckets = []
                self.Buckets = make_buckets(self.BottomLeft, self.TopRight, 
                                            self.AllBuckets)

                for el, el_bbox in self.Elements:
                    self.insertIntoSubdivision(el, el_bbox)

                # Free up some memory
                del self.Elements

                self.insertIntoSubdivision(element, bbox)
            else:
                # Simple:
                self.Elements.append((element, bbox))
        else:
            self.insertIntoSubdivision(element, bbox)

    def insertIntoSubdivision(self, element, bbox):
        # Surprisingly, this can't be done much faster.
        for bucket in self.AllBuckets:
            if do_boxes_intersect((bucket.bottom_left(), bucket.top_right()), bbox):
                bucket.add_element(element, bbox)

    def find_element(self, point):
        if self.Buckets:
            # We have subdivisions. Use them.
            (dimensions,) = point.shape
            bucket = self.Buckets
            for dim in range(0, dimensions):
                if point[ dim ] < self.Center[ dim ]:
                    bucket = bucket[ 0 ]
                else:
                    bucket = bucket[ 1 ]

            return bucket.find_element(point)
        else:
            # We don't. Perform linear search.
            for el, bbox in self.Elements:
                if el.contains_point(point):
                    return el
            raise RuntimeError, "No element found at %s." % str(point)

    def visualize(self, file):
        file.write("%f %f\n" % (self.BottomLeft[0], self.BottomLeft[1]));
        file.write("%f %f\n" % (self.TopRight[0], self.BottomLeft[1]));
        file.write("%f %f\n" % (self.TopRight[0], self.TopRight[1]));
        file.write("%f %f\n" % (self.BottomLeft[0], self.TopRight[1]));
        file.write("%f %f\n\n" % (self.BottomLeft[0], self.BottomLeft[1]));
        if self.Buckets:
            for i in self.AllBuckets:
                i.visualize(file)




def build_spatial_binary_tree(elements):
    elements_with_boxes = map(lambda el: (el, el.bounding_box()), elements)
    bottom_left,top_right = _get_elements_bounding_box(elements_with_boxes)
    bucket = SpatialBinaryTreeBucket(bottom_left, top_right)
    for el, bbox in elements_with_boxes:
        bucket.add_element(el, bbox)

        # test
        if False:
            center = pytools.average([node.Coordinates for node in el.nodes()])
            assert el.contains_point(center)
            bucket.find_element(center)
    return bucket




def build_element_finder(elements):
    bucket = build_spatial_binary_tree(elements)
    return bucket.find_element

