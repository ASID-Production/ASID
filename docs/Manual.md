<h1 align="center">ASID Manual</h1>

1. [Points menu](#Points_menu)
   1. Points
   2. Points lists
   3. Context menu
      * Attach representation
      * Detach Representation
      * Add list
      * Add point
      * Delete
   4. Reordering
   5. Selection
2. [Properties menu](#Properties_menu)
   1. Editing
      * Adding new fields
      * Deleting fields
      * Multiple editing
3. [Visual space](#Visual_space)
   * Controls
     * Rotation
     * Translation
     * Selection
4. Uniforms menu
5. Representations
6. Extensions
   1. Find sub
   2. DB Search
   3. Open
<a name="Points_menu"></a>
### 1. Points menu
#### 1.1 Points
Minimal structural element, can hold data and will be used by representations.
#### 1.2 Points lists
Container of points, can hold data, but will not be used by representations.
#### 1.3 Context menu
##### Attach representation
In this menu you can attach representation to point or point list, attaching representation to list is similar to attaching representation to every point in list individually. Points or list with representation will be visible in visual space.
##### Detach representation
In this menu you can detach representation from point or point list if there is any, when detached point will no longer be visible in visual space.
##### Add list
Create new empty list of points.
##### Add point
Create new point at parent list.
##### Delete
Delete a point or list with it content. All dependencies will be replaced with the values set at the time of deletion
#### 1.4 Reordering
Some representations are order sensitive. In that case reordering can be preformed by drag'n'drop of a point in within the native list. 
#### 1.5 Selection
Multiple lists and points can be selected by holding shift or ctrl button for multiple editing.
<a name="Properties_menu"></a>
### 2. Properties menu
Properties menu contain info about attached to point or point list data (fields).
When multiple points/lists selected all equal fields will be shown as usual, data of different fields will be shown as blank line.
Properties can be inherited by special syntax @seq (example: @123, property will be inherited from point/list with seq property == 123)
'None' is reserved value and mean no value, signal for using default values in representations.
Seq property can not be redacted, its unique id of point/list.
Some properties can be used by representations for visualization, see [Representations]() for more info.
#### 2.1 Editing
Double tap on a property to edit it value. You can not edit property name, this will be ignored, and you don't need to specify its name.
Special syntax for list of floats - [f,f,f] (list of 3 floats), inheritance - @seq, reserved word 'None', without quotes.
##### Adding new properties
When selected, type name of a property, then it will be created with default value 'None'.
##### Deleting properties
Property will be removed from point/list.
##### Multiple editing
You can set property value for multiple points/lists at once.
<a name="Visual_Space"></a>
### 3. Visual Space
#### Controls
##### Rotation
Left or right mouse buttons
##### Translation
Shift+left mouse button
##### Selection
Select - left mouse button, deselect - right mouse button
### Uniforms menu
Variables of visual space, can be edited like properties
### Representations
property - format (default value)
#### Sphere
Sphere representation, searching for properties:
* coord - [f,f,f] ([0,0,0])
* color - [f,f,f,f] ([0,0,0,1])
* rad - f (1)
* pick - f (0)
#### Bond
Bond representation (hole tube, order sensitive, uses two points in order), searching for properties:
* coord - [f,f,f] ([0,0,0])
* color - [f,f,f,f] ([0,0,0,1])
* rad - f (1)
* pick - f (0)
#### Label
Label representation, searching for properties:
* coord - [f,f,f] ([0,0,0])
* label - string (None)
#### Plane
Plane representation (single polygon, order sensitive, uses three points in order), searching for properties:
* coord - [f,f,f] ([0,0,0])
* color - [f,f,f,f] ([0,0,0,1])
#### Line
Line representation (order sensitive, uses two points in order), searching for properties:
* coord - [f,f,f] ([0,0,0])
* color - [f,f,f,f] ([0,0,0,1])
* pick - f (0)
#### Dashed line
Dashed line representation (order sensitive, uses two points in order), searching for properties:
* coord - [f,f,f] ([0,0,0])
* color - [f,f,f,f] ([0,0,0,1])
* pick - f (0)
