import xml.etree.ElementTree as ET
import os
from PIL import Image, ImageOps
import hashlib
import argparse

class TileOp:
    def __init__(self, tile: "Tile", flip_h: bool = False, flip_v: bool = False, transpose: bool = False) -> None:
        self.tile = tile
        self.flip_h = flip_h
        self.flip_v = flip_v
        self.transpose = transpose

    def __repr__(self):
        transformations = []
        if self.flip_h:
            transformations.append("flip_h")
        if self.flip_v:
            transformations.append("flip_v")
        if self.transpose:
            transformations.append("transpose")
        
        transformations_str = ", ".join(transformations) if transformations else "None"
        
        return f"TileOp(tile=Tile(ID={self.tile.tile_id}), transformations=[{transformations_str}])"

class Tile:
    def __init__(self, tile_id, image, x_ref, y_ref):
        self.tile_id = tile_id  # Unique identifier of the tile
        self.image = image  # Image of the tile
        self.x_ref = x_ref  # x position of the tile in the assembled image
        self.y_ref = y_ref  # y position of the tile in the assembled image

    def __repr__(self):
        return f"Tile(ID={self.tile_id}, x_ref={self.x_ref}, y_ref={self.y_ref})"
    
    def generate_transformations(self):
        """
        Generates possible transformations for a tile:
        - flips: horizontal, vertical, and diagonal (swapping x and y axes)
        Returns a list of tuples (transformed image, flip horizontal, flip vertical, flip diagonal)
        """
        transformations = [
            (self.image, TileOp(self, False, False, False)),  # No flip
            (ImageOps.mirror(self.image), TileOp(self, True, False, False)),  # Horizontal flip
            (ImageOps.flip(self.image), TileOp(self, False, True, False)),  # Vertical flip
            (self.image.transpose(Image.Transpose.TRANSPOSE), TileOp(self, False, False, True)),  # Diagonal flip (swap x and y)
            (ImageOps.mirror(self.image.transpose(Image.Transpose.TRANSPOSE)), TileOp(self, True, False, True)),  # Horizontal + diagonal flip
            (ImageOps.flip(self.image.transpose(Image.Transpose.TRANSPOSE)), TileOp(self, False, True, True)),  # Vertical + diagonal flip
        ]
        return transformations


class Cell:
    def __init__(self, tile_op: TileOp, x, y):
        self.tile_op = tile_op  # Reference to an instance of the Tile class
        self.x = x  # x position of the cell in the original image
        self.y = y  # y position of the cell in the original image

    def __repr__(self):
        return (f"Cell(TileOp={self.tile_op}, x={self.x}, y={self.y})")

def hash_image(image):
    """
    Creates a unique hash for an image using hashlib.
    """
    image_data = image.tobytes()
    return hashlib.md5(image_data).hexdigest()

def split_image_into_tiles(image_path, tile_size):
    # Load the image
    image = Image.open(image_path)
    width, height = image.size

    # Calculate the number of tiles along the x and y axes
    num_tiles_x = width // tile_size
    num_tiles_y = height // tile_size

    # List to store unique tiles (and their transformations)
    unique_tiles = []
    
    # Dictionary of transformations (key: hash, value: reference to a Tile)
    tile_hashes = {}

    # List to represent the final image with transformations
    image_representation = []
    total_tiles = num_tiles_x * num_tiles_y
    tile_count = 0
    # Split the image into tiles
    for i in range(num_tiles_y):
        for j in range(num_tiles_x):

            # Calculate the coordinates of each tile
            left = j * tile_size
            upper = i * tile_size
            right = left + tile_size
            lower = upper + tile_size

            # Crop the tile
            tile_image = image.crop((left, upper, right, lower))

            # Calculate the hash of the original tile
            h = hash_image(tile_image)

            # Check if any of the transformations of this tile already exist
            if h in tile_hashes:
                # If the hash exists, retrieve the corresponding tile
                tile_op = tile_hashes[h]
            else:
                # Otherwise, generate all transformations and add their hashes to the dictionary
                tile_id = len(unique_tiles) + 1
                new_tile = Tile(tile_id=tile_id, image=tile_image, x_ref=j, y_ref=i)
                unique_tiles.append(new_tile)
                transformations = new_tile.generate_transformations()
                # Save all transformations in the dictionary
                for transf_image, op in transformations:
                    h_transf = hash_image(transf_image)
                    tile_hashes[h_transf] = op
                tile_op = transformations[0][1]
                # Add the cell with the original tile
            cell_rep = Cell(tile_op=tile_op, x=j, y=i)
            image_representation.append(cell_rep)
            tile_count += 1
            if tile_count % 1000 == 0:
                print(f"Tile {tile_count} out of {total_tiles}, unique {len(unique_tiles)}")

    return image_representation, unique_tiles

def assemble_tiles(unique_tiles, tile_size):
    """
    Assemble all unique tiles into a large square image.
    """
    # Calculate the grid size to display the tiles
    grid_size = int(len(unique_tiles) ** 0.5) + 1

    # Create a large empty square image to host all tiles
    large_image = Image.new('RGBA', (grid_size * tile_size, grid_size * tile_size))
    x = 0
    y = 0
    # Place each unique tile in the large image
    for tile in unique_tiles:
        x_pos = x * tile_size
        y_pos = y * tile_size
        x += 1
        if x > grid_size - 1:
            x = 0
            y += 1
        position = (x_pos, y_pos)
        large_image.paste(tile.image, position)

    return large_image


def encode_tile(tile_id, flip_h, flip_v, flip_diag):
    """
    Encode a tile with its ID, horizontal flip, vertical flip, and diagonal flip.
    Flips are encoded in bits 29, 30, and 31.
    """
    # Start with the tile ID
    encoded_value = tile_id
    
    # Horizontal flip is encoded in bit 31
    if flip_h:
        encoded_value |= (1 << 31)  # Bit 31 for horizontal flip
    
    # Vertical flip is encoded in bit 30
    if flip_v:
        encoded_value |= (1 << 30)  # Bit 30 for vertical flip
    
    # Diagonal flip (90° rotation or swapping x and y axes) is encoded in bit 29
    if flip_diag:
        encoded_value |= (1 << 29)  # Bit 29 for diagonal flip
    
    return encoded_value


def generate_tiled_xml(unique_tiles, image_representation, tile_size, map_width, map_height, output_file="map.tmx", image_tileset_name="", layer_name="Graphics"):
    """
    Generates an XML file compatible with Tiled using the given tiles and image representation.
    """
    # Create the root element of the map
    map_elem = ET.Element('map', {
        'version': '1.10',
        'tiledversion': '1.11.0',
        'orientation': 'orthogonal',
        'renderorder': 'left-up',
        'width': str(map_width),
        'height': str(map_height),
        'tilewidth': str(tile_size),
        'tileheight': str(tile_size),
        'infinite': '0',
        'nextlayerid': '2',
        'nextobjectid': '1'
    })
    
    cols = int(len(unique_tiles) ** 0.5) + 1

    # Add the tileset element
    tileset_elem = ET.SubElement(map_elem, 'tileset', {
        'firstgid': '1',
        'name': 'assembled_tiles',
        'tilewidth': str(tile_size),
        'tileheight': str(tile_size),
        'tilecount': str(cols**2),
        'columns': str(cols)
    })
    # Image source for the tileset
    ET.SubElement(tileset_elem, 'image', {
        'source': image_tileset_name,
        'width': str(int(cols) * tile_size),
        'height': str(int(cols) * tile_size)
    })

    # Add a layer for the tiles
    layer_elem = ET.SubElement(map_elem, 'layer', {
        'id': '1',
        'name': layer_name,
        'width': str(map_width),
        'height': str(map_height)
    })

    # Add the layer data in CSV format
    data_elem = ET.SubElement(layer_elem, 'data', {'encoding': 'csv'})

    # Generate the encoded data for each tile and add it as CSV
    csv_data = []
    for c in image_representation:
        tile_id = c.tile_op.tile.tile_id
        encoded_value = encode_tile(tile_id, c.tile_op.flip_h, c.tile_op.flip_v, c.tile_op.transpose)
        csv_data.append(str(encoded_value))
    
    # Format as rows to match the map width
    rows = []
    for i in range(0, len(csv_data), map_width):
        rows.append(",".join(csv_data[i:i + map_width]))
    data_elem.text = ",\n".join(rows)

    # Write the XML file
    tree = ET.ElementTree(map_elem)
    with open(output_file, "wb") as files:
        tree.write(files, encoding="UTF-8", xml_declaration=True)


def image_to_tiled(image_path: str, tile_size: int = 8):
    image_representation, unique_tiles = split_image_into_tiles(image_path, tile_size)

    # Assemble all unique tiles into a large square image
    large_image = assemble_tiles(unique_tiles, tile_size)

    fname = "joined_tiles_" + image_path

    # Save the large image containing all unique tiles
    large_image.save(fname)

    last_cell = image_representation[-1]
    print(f"There are {len(unique_tiles)} unique tiles.")
    
    # Generate the Tiled file
    map_width = last_cell.x + 1  # Example map width
    map_height = last_cell.y + 1  # Example map height
    file_name, _ = os.path.splitext(image_path)
    map_name = file_name + ".tmx"
    generate_tiled_xml(unique_tiles, image_representation, tile_size, map_width, map_height, output_file=map_name, image_tileset_name=fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert an image into Tiled format.')
    parser.add_argument('image_path', type=str, help='Path to the input image')
    parser.add_argument('--tile_size', type=int, default=8, help='Size of the tiles (default: 8)')
    args = parser.parse_args()

    image_to_tiled(args.image_path, args.tile_size)
