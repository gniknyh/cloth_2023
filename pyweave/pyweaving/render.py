from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from math import asin, cos, floor, sin

import os.path

from PIL import Image, ImageDraw, ImageFont

from pyweaving import WarpThread, WeftThread

import numpy as np
import imageio
import os
os.environ["OPENCV_IO_ENABLE_OPENEXR"]="1"
import cv2
psi = 0.523
psi = 0.0
# u_max = 0.523 ### pi/6  30 deg
# u_max = 0.0314 ### pi/6  30 deg
u_max = 0.314 ### pi/6  30 deg
# u_max = 0.5 * 2
sin_u_max = sin(u_max)

sin_psi = sin(psi)
cos_psi = cos(psi)

__here__ = os.path.dirname(__file__)

font_path = os.path.join(__here__, 'data', 'Arial.ttf')


class ImageRenderer(object):
    # TODO:
    # - Add a "drawndown only" option
    # - Add a default tag (like a small delta symbol) to signal the initial
    # shuttle direction
    # - Add option to render the backside of the fabric
    # - Add option to render a bar graph of the thread crossings along the
    # sides
    # - Add option to render 'stats table'
    #   - Number of warp threads
    #   - Number of weft threads
    #   - Number of harnesses/shafts
    #   - Number of treadles
    #   - Warp unit size / reps
    #   - Weft unit size / reps
    #   - Longest warp float
    #   - Longest weft float
    #   - Selvedge continuity
    # - Add option to rotate orientation
    # - Add option to render selvedge continuity
    # - Add option to render inset "scale view" rendering of fabric
    # - Add option to change thread spacing
    # - Support variable thickness threads
    # - Add option to render heddle count on each shaft
    def __init__(self, draft, liftplan=None, margin_pixels=20, scale=20,
                 foreground=(127, 127, 127), background=(255, 255, 255),
                 markers=(0, 0, 0), numbering=(200, 0, 0)):
        self.draft = draft

        self.liftplan = liftplan

        self.margin_pixels = margin_pixels
        self.pixels_per_square = scale

        self.background = background
        self.foreground = foreground
        self.markers = markers
        self.numbering = numbering

        self.font_size = int(round(scale * 1.2))

        self.font = ImageFont.truetype(font_path, self.font_size)

    def pad_image(self, im):
        w, h = im.size
        desired_w = w + (self.margin_pixels * 2)
        desired_h = h + (self.margin_pixels * 2)
        new = Image.new('RGB', (desired_w, desired_h), self.background)
        new.paste(im, (self.margin_pixels, self.margin_pixels))
        return new

    def make_pil_image(self):
        # width_squares = len(self.draft.warp) + 6
        # if self.liftplan or self.draft.liftplan:
        #     width_squares += len(self.draft.shafts)
        # else:
        #     width_squares += len(self.draft.treadles)

        # height_squares = len(self.draft.weft) + 6 + len(self.draft.shafts)
        width_squares = len(self.draft.warp)
        height_squares = len(self.draft.weft)

        # XXX Not totally sure why the +1 is needed here, but otherwise the
        # contents overflows the canvas
        width = (width_squares * self.pixels_per_square)
        height = (height_squares * self.pixels_per_square)

        im = Image.new('RGB', (width, height), self.background)
        self.normal_map = Image.new('RGB', (width, height), self.background)
        self.tangent_map = Image.new('RGB', (width, height), self.background)

        draw = ImageDraw.Draw(im) 

        # self.paint_warp(draw)
        # self.paint_threading(draw)

        # self.paint_weft(draw)
        # if self.liftplan or self.draft.liftplan:
        #     self.paint_liftplan(draw)
        # else:
        #     self.paint_tieup(draw)
        #     self.paint_treadling(draw)
        # self.paint_tieup(draw)
        # self.paint_treadling(draw)
        self.paint_drawdown(draw, im )
        # self.paint_start_indicator(draw)
        del draw 

        # im = self.pad_image(im)
        return im 

    # def paint_start_indicator(self, draw):
    #     endy = ((len(self.draft.shafts) + 6) * self.pixels_per_square) - 1
    #     starty = (endy - (self.pixels_per_square // 2))
    #     if self.draft.start_at_lowest_thread:
    #         # right side
    #         endx = len(self.draft.warp) * self.pixels_per_square
    #         startx = endx - self.pixels_per_square
    #     else:
    #         # left side
    #         startx = 0
    #         endx = self.pixels_per_square
    #     vertices = [
    #         (startx, starty),
    #         (endx, starty),
    #         (startx + (self.pixels_per_square // 2), endy),
    #     ]
    #     draw.polygon(vertices, fill=self.markers)

    def paint_warp(self, draw):
        starty = 0
        endy = self.pixels_per_square
        for ii, thread in enumerate(self.draft.warp):
            # paint box, outlined with foreground color, filled with thread
            # color
            startx = self.pixels_per_square * ii
            endx = startx + self.pixels_per_square
            draw.rectangle((startx, starty, endx, endy),
                           outline=self.foreground,
                           fill=thread.color.rgb)

    def paint_fill_marker(self, draw, box):
        startx, starty, endx, endy = box
        draw.rectangle((startx + 2, starty + 2, endx - 2, endy - 2),
                       fill=self.markers)

    def paint_threading(self, draw):
        num_threads = len(self.draft.warp)
        num_shafts = len(self.draft.shafts)

        for ii, thread in enumerate(self.draft.warp):
            startx = (num_threads - ii - 1) * self.pixels_per_square
            endx = startx + self.pixels_per_square

            for jj, shaft in enumerate(self.draft.shafts):
                starty = (4 + (num_shafts - jj)) * self.pixels_per_square
                endy = starty + self.pixels_per_square
                draw.rectangle((startx, starty, endx, endy),
                               outline=self.foreground)

                if shaft == thread.shaft:
                    # draw threading marker
                    self.paint_fill_marker(draw, (startx, starty, endx, endy))

            # paint the number if it's a multiple of 4
            thread_no = ii + 1
            if ((thread_no != num_threads) and
                (thread_no != 0) and
                    (thread_no % 4 == 0)):
                # draw line
                startx = endx = (num_threads - ii - 1) * self.pixels_per_square
                starty = 3 * self.pixels_per_square
                endy = (5 * self.pixels_per_square) - 1
                draw.line((startx, starty, endx, endy),
                          fill=self.numbering)
                # draw text
                draw.text((startx + 2, starty + 2),
                          str(thread_no),
                          font=self.font,
                          fill=self.numbering)

    def paint_weft(self, draw):
        offsety = (6 + len(self.draft.shafts)) * self.pixels_per_square
        startx_squares = len(self.draft.warp) + 5
        if self.liftplan or self.draft.liftplan:
            startx_squares += len(self.draft.shafts)
        else:
            startx_squares += len(self.draft.treadles)
        startx = startx_squares * self.pixels_per_square
        endx = startx + self.pixels_per_square

        for ii, thread in enumerate(self.draft.weft):
            # paint box, outlined with foreground color, filled with thread
            # color
            starty = (self.pixels_per_square * ii) + offsety
            endy = starty + self.pixels_per_square
            draw.rectangle((startx, starty, endx, endy),
                           outline=self.foreground,
                           fill=thread.color.rgb)

    # def paint_liftplan(self, draw):
    #     num_threads = len(self.draft.weft)

    #     offsetx = (1 + len(self.draft.warp)) * self.pixels_per_square
    #     offsety = (6 + len(self.draft.shafts)) * self.pixels_per_square
    #     for ii, thread in enumerate(self.draft.weft):
    #         starty = (ii * self.pixels_per_square) + offsety
    #         endy = starty + self.pixels_per_square

    #         for jj, shaft in enumerate(self.draft.shafts):
    #             startx = (jj * self.pixels_per_square) + offsetx
    #             endx = startx + self.pixels_per_square
    #             draw.rectangle((startx, starty, endx, endy),
    #                            outline=self.foreground)

    #             if shaft in thread.connected_shafts:
    #                 # draw liftplan marker
    #                 self.paint_fill_marker(draw, (startx, starty, endx, endy))

    #         # paint the number if it's a multiple of 4
    #         thread_no = ii + 1
    #         if ((thread_no != num_threads) and
    #             (thread_no != 0) and
    #                 (thread_no % 4 == 0)):
    #             # draw line
    #             startx = endx
    #             starty = endy
    #             endx = startx + (2 * self.pixels_per_square)
    #             endy = starty
    #             draw.line((startx, starty, endx, endy),
    #                       fill=self.numbering)
                          
    def paint_tieup(self, draw):
        offsetx = (1 + len(self.draft.warp)) * self.pixels_per_square
        offsety = 5 * self.pixels_per_square

        num_treadles = len(self.draft.treadles)
        num_shafts = len(self.draft.shafts)

        for ii, treadle in enumerate(self.draft.treadles):
            startx = (ii * self.pixels_per_square) + offsetx
            endx = startx + self.pixels_per_square

            treadle_no = ii + 1

            for jj, shaft in enumerate(self.draft.shafts):
                starty = (((num_shafts - jj - 1) * self.pixels_per_square) +
                          offsety)
                endy = starty + self.pixels_per_square

                draw.rectangle((startx, starty, endx, endy),
                               outline=self.foreground)

                if shaft in treadle.shafts:
                    self.paint_fill_marker(draw, (startx, starty, endx, endy))

                # on the last treadle, paint the shaft markers
                if treadle_no == num_treadles:
                    shaft_no = jj + 1
                    if (shaft_no != 0) and (shaft_no % 4 == 0):
                        # draw line
                        line_startx = endx
                        line_endx = line_startx + (2 * self.pixels_per_square)
                        line_starty = line_endy = starty
                        draw.line((line_startx, line_starty,
                                   line_endx, line_endy),
                                  fill=self.numbering)
                        draw.text((line_startx + 2, line_starty + 2),
                                  str(shaft_no),
                                  font=self.font,
                                  fill=self.numbering)

            # paint the number if it's a multiple of 4 and not the first one
            if (treadle_no != 0) and (treadle_no % 4 == 0):
                # draw line
                startx = endx = (treadle_no * self.pixels_per_square) + offsetx
                starty = 3 * self.pixels_per_square
                endy = (5 * self.pixels_per_square) - 1
                draw.line((startx, starty, endx, endy),
                          fill=self.numbering)
                # draw text on left side, right justified
                textw, texth = draw.textsize(str(treadle_no), font=self.font)
                draw.text((startx - textw - 2, starty + 2),
                          str(treadle_no),
                          font=self.font,
                          fill=self.numbering)

    def paint_treadling(self, draw):
        num_threads = len(self.draft.weft)

        offsetx = (1 + len(self.draft.warp)) * self.pixels_per_square
        offsety = (6 + len(self.draft.shafts)) * self.pixels_per_square
        for ii, thread in enumerate(self.draft.weft):
            starty = (ii * self.pixels_per_square) + offsety
            endy = starty + self.pixels_per_square

            for jj, treadle in enumerate(self.draft.treadles):
                startx = (jj * self.pixels_per_square) + offsetx
                endx = startx + self.pixels_per_square
                draw.rectangle((startx, starty, endx, endy),
                               outline=self.foreground)

                if treadle in thread.treadles:
                    # draw treadling marker
                    self.paint_fill_marker(draw, (startx, starty, endx, endy))

            # paint the number if it's a multiple of 4
            thread_no = ii + 1
            if ((thread_no != num_threads) and
                (thread_no != 0) and
                    (thread_no % 4 == 0)):
                # draw line
                startx = endx
                starty = endy
                endx = startx + (2 * self.pixels_per_square)
                endy = starty
                draw.line((startx, starty, endx, endy),
                          fill=self.numbering)
                # draw text
                draw.text((startx + 2, starty - 2 - self.font_size),
                          str(thread_no),
                          font=self.font,
                          fill=self.numbering)

    def normal_to_rgb_int(self, normal):
        pix_normal_r = floor(127.5*(normal[0]+1))
        pix_normal_g = floor(127.5*(normal[1]+1))
        pix_normal_b = floor(127.5*(normal[2]+1))
        pix_normal_rgb = (pix_normal_r, pix_normal_g, pix_normal_b)
        return pix_normal_rgb

    def normal_to_rgb_gamma_correct(self, normal):
        pix_normal_r = ((normal[0]+1)*0.5)** ( 2.2)
        pix_normal_g = ((normal[1]+1)*0.5)** ( 2.2)
        pix_normal_b = ((normal[2]+1)*0.5)** ( 2.2) 
        pix_normal_rgb = (pix_normal_r, pix_normal_g, pix_normal_b)
        return pix_normal_rgb
 

    def paint_drawdown(self, draw, im):
        offsety = 0
        floats = self.draft.compute_floats()

        img_size = (self.normal_map.size[1],self.normal_map.size[0], 3)
        normal_mat = np.zeros(img_size)
        tangent_mat = np.zeros(img_size)

        normal_mat_float = np.zeros(img_size)
        tangent_mat_float = np.zeros(img_size)

        for start, end, visible, length, thread in floats:
            if visible:
                startx = start[0] * self.pixels_per_square
                starty = (start[1] * self.pixels_per_square) + offsety
                endx = (end[0] + 1) * self.pixels_per_square
                endy = ((end[1] + 1) * self.pixels_per_square) + offsety

                draw.rectangle((startx, starty, endx, endy),
                               outline=self.foreground,
                               fill=thread.color.rgb)

                center_x = (startx + endx)/2
                center_y = (starty + endy)/2                

                for x in range(startx, endx):
                    for y in range(starty, endy):
                        pix_x = x + 0.5
                        pix_y = y + 0.5
                        ratio_y_l = 2*(center_y-pix_y)/(endy-starty)
                        ratio_x_w = 2*(pix_x-center_x)/(endx-startx)
                        
                        u = asin(ratio_y_l * sin_u_max)
                        v = asin(ratio_x_w)
                        if isinstance(thread, WeftThread):
                            v = asin(ratio_x_w * sin_u_max)
                            u = asin(ratio_y_l)

                        sin_u = sin(u); sin_v = sin(v); cos_u = cos(u); cos_v = cos(v)
                        pix_normal = (sin_v, sin_u*cos_v, cos_u*cos_v)
                        pix_tangent = (
                            -cos_v*sin_psi, 
                            cos_u*cos_psi + sin_u*sin_v*sin_psi,
                            -sin_u*cos_psi + cos_u*sin_v*sin_psi

                        )
                        
                        pix_normal_rgb = self.normal_to_rgb_int(pix_normal)
                        pix_tangent_rgb = self.normal_to_rgb_int(pix_tangent) 

                        pix_normal_rgb_f = self.normal_to_rgb_gamma_correct(pix_normal)
                        pix_tangent_rgb_f = self.normal_to_rgb_gamma_correct(pix_tangent) 
                        
                        normal_mat[y,x,2] = pix_normal_rgb[0]
                        normal_mat[y,x,1] = pix_normal_rgb[1]
                        normal_mat[y,x,0] = pix_normal_rgb[2]
                        tangent_mat[y,x,2] = pix_tangent_rgb[0]
                        tangent_mat[y,x,1] = pix_tangent_rgb[1]
                        tangent_mat[y,x,0] = pix_tangent_rgb[2]

                        normal_mat_float[y,x,2] = pix_normal_rgb_f[0]
                        normal_mat_float[y,x,1] = pix_normal_rgb_f[1]
                        normal_mat_float[y,x,0] = pix_normal_rgb_f[2]
                        tangent_mat_float[y,x,2] = pix_tangent_rgb_f[0]
                        tangent_mat_float[y,x,1] = pix_tangent_rgb_f[1]
                        tangent_mat_float[y,x,0] = pix_tangent_rgb_f[2]
                        
        
        cv2.imwrite('normal.png', normal_mat.astype(np.float32))
        cv2.imwrite('tangent_map.png', tangent_mat.astype(np.float32))

        cv2.imwrite('normal.exr', normal_mat_float.astype(np.float32))
        cv2.imwrite('tangent_map.exr', tangent_mat_float.astype(np.float32))
        
    def set_pixel_resolution(self, value):
        self.pixels_per_square = value

    def show(self):
        im  = self.make_pil_image()
        im.show()

    def save(self, filename):
        im = self.make_pil_image()
        im.save(filename)
        


