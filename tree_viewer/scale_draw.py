from PyQt4 import QtGui
from ete2.treeview.qt4_render import _EmptyItem


def add_scale(img, mainRect, parent, length, scale_type, min_value, max_value):
    # img.scale_length, img.scale_type, img.scale_min_value, img.scale_max_value
    if img.show_scale:
        scale_item = _EmptyItem()
        custom_pen = QtGui.QPen(QtGui.QColor("black"), 1)

        line = QtGui.QGraphicsLineItem(scale_item)
        line.setPen(custom_pen)
        calculate_marks = length/10
        mark = {}
        for x in range(calculate_marks):
            mark[x] = QtGui.QGraphicsLineItem(scale_item)
            mark[x].setPen(custom_pen)

        line.setLine(0, 5, length, 5)
        for x in range(calculate_marks):
            mark[x].setLine(0, 0, 0, 10)


        # line3.setLine(length, 0, length, 10)

        length_text = float(length) / img._scale if img._scale else 0.0
        scale_text = "%0.2f" % (length_text)
        scale = QtGui.QGraphicsSimpleTextItem(scale_text)
        scale.setParentItem(scale_item)
        scale.setPos(0, 10)

        scale_item.setParentItem(parent)
        dw = max(0, length-mainRect.width())
        scale_item.setPos(mainRect.bottomLeft())
        scale_item.moveBy(img.margin_left, 0)
        mainRect.adjust(0, 0, dw, length)
